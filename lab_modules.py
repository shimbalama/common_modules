#!/usr/bin/env python3

import sys
from subprocess import call, Popen, PIPE
import os
from collections import defaultdict
from Bio import SeqIO
from glob import glob
import shutil

def main():
    pass

def call_iqtree_with_SNP_aln(args, d_count):

    constants=','.join([str(d_count.get('A')), str(d_count.get('C')),
           str(d_count.get('G')), str(d_count.get('T'))])
    print ('constants',constants)
    call(['nice', 'iqtree', '-s', args.out + '_trimmed_SNP.aln',
          '-m', 'GTR+G', '-bb', '1000', '-czb', '-fconst', constants])


def get_seq_type(seqs):

    #better with regex?
    amino_acids = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    nucleotides = ['A','T','C','G','N','R','Y','S','W','K','M','B','D','H','V','N']
    # from http://www.bioinformatics.org/sms/iupac.html
    protein = False
    DNA = True
    b = False
    for record in SeqIO.parse(seqs,'fasta'):
        for pos in str(record.seq).upper():
            if pos not in nucleotides:
                DNA  = False
                protein = True
                b = True
                break #speed it up..
            elif pos not in nucleotides and pos not in amino_acids:
                print ('Error, not DNA or protein: ', pos, ' in ', record)
                protein = False
                sys.exit(0)
        break
    if DNA:
        seq_type = 'DNA'
    elif protein:
        seq_type = 'prot'
    return seq_type

def cat(args):

    '''
    Concatenates assemblies; takes n assemblies and makes one file with all contigs
    '''
    print ('cating assemblies...')
    if os.path.isfile(args.input_folder + '/concatenated_assemblies/concatenated_assemblies.fasta'):
        make = False
    else:
        if not os.path.exists(args.input_folder + '/concatenated_assemblies'):
            os.makedirs(args.input_folder + '/concatenated_assemblies')
        make = True
        fout = open(args.input_folder + '/concatenated_assemblies/concatenated_assemblies.fasta', 'w')
    assemblies = set([])
    if make:
        for assembly in glob(args.input_folder + '/*'):
            ass_file_name = assembly.strip().split('/')[-1].split('.')
            if ass_file_name == ['concatenated_assemblies']:
                continue
            try:
                assert len(ass_file_name) == 2
            except:
                print ('Please ensure that your database file name only has one fullstop i.e., sample_1.fa NOT sample.1.fa')
                sys.exit(0)
            ass_file_name = ass_file_name[0].replace('#','_')#Hashes break stuff
            assemblies.add(ass_file_name)
            for i, record in enumerate(SeqIO.parse(assembly, 'fasta')):
                record.id = 'gnl|MYDB|'+ass_file_name + '_' + str(i)#can't handle all the variablity any other way
                SeqIO.write(record,fout,'fasta')
    if make:
        fout.close()
    #If pre cated
    print ('Getting assembly names...')
    for record in SeqIO.parse(args.input_folder + '/concatenated_assemblies/concatenated_assemblies.fasta', 'fasta'):
        ass = '_'.join(record.id.split('_')[:-1]).replace('gnl|MYDB|','')#just takes contig number off the end
        if make:
            assert ass in assemblies
        else:
            assemblies.add(ass)
    print (len(assemblies), 'assemblies used')

    return assemblies

def fasta(tup):
    
    '''
    Makes a fasta file of all hits for each query
    '''
    args, query = tup 
    #make oonfig file from blast output
    current = query + '_config.txt'
    fout = open(current,'w')
    coords = []
    with open('blast_out.txt', 'r') as fin:
        for line in fin:
            bits = line.strip().split()
            if float(bits[2]) < (float(args.percent_identity) -1.0):
                continue
            if bits[0] == query:
                contig = bits[1]
                seq_start = bits[8]
                seq_end = bits[9]
                coords.append(seq_start + '-' + seq_end) #only want in one direction as just an id 
                if int(seq_start) > int(seq_end):#minus/plus is right, although seems backwards..
                    line = ' '.join([contig, seq_end + '-' + seq_start, 'minus'])
                    fout.write(line + '\n')
                else:
                    line = ' '.join([contig, seq_start + '-' + seq_end, 'plus'])
                    fout.write(line + '\n')
    fout.close()
    #Cut seq of interest out of contig
    cat = args.input_folder + '/concatenated_assemblies/concatenated_assemblies.fasta'
    cmd = ['blastdbcmd',
        '-db', cat,
        '-entry_batch', current,
        '-out', query + '_all_nuc_seqs_tmp.fasta']
    call(cmd)
    
    #add coords to fasta as unique id for multiple hits in same contig
    with open(query + '_all_nuc_seqs.fasta','w') as fout:
        for i, record in enumerate(SeqIO.parse(query + '_all_nuc_seqs_tmp.fasta','fasta')):
            record.description = coords[i]
            record.id = str(record.id).split(':')[1]
            seq = str(record.seq)
            seq_len = len(seq)
            if seq_len%3==0:
                SeqIO.write(record,fout,'fasta')          
            else:
                while seq_len%3!=0:
                    seq_len -= 1
                seq = seq[:seq_len]
                try:assert len(seq)%3==0
                except: print (len(seq), len(seq)%3)
                fout.write('>'+record.id + ' ' + record.description +'\n')
                for i in range(0, seq_len, 60):
                    fout.write(seq[i:i+60] +'\n')            
    

def get_query_seqs(args):

    query_seqs = {}
    for record in SeqIO.parse(args.query, 'fasta'):
        query_seqs[record.id] = record

    return query_seqs

def blast(args):

    '''
    Blasts (Nuc or Prot) seq/s of interest (Query) against db of concatenated
    assemblies(contigs)/CDS/proteins
    '''
    print ('testing query seqs')
    seq_type_query = get_seq_type(args.query)
    print ('fine -> ', seq_type_query)
    print ('testing db seqs')
    cat = args.input_folder + '/concatenated_assemblies/concatenated_assemblies.fasta'
    seq_type_db = get_seq_type(cat)
    print ('fine -> ', seq_type_db)
    #db
    if seq_type_db == 'DNA':
        seq_type = 'nucl'
    else:
        seq_type = 'prot'
    try:
        call(['makeblastdb',
        '-in', cat,
        '-parse_seqids',
        '-dbtype', seq_type])#'-parse_seqids',
    except:
        print ('you need to install makeblastdb or fix paths to it')

    #blast
    if seq_type_query == 'DNA' and seq_type_db == 'DNA':
        blast_type = 'blastn'
        call_blast(args, cat, 'blastn')
        print ('Running blastn...')
    elif seq_type_query == 'prot' and seq_type_db == 'prot':
        blast_type = 'blastp'
        call_blast(args, cat, 'blastp')
        print ('Running blastp...')
    elif seq_type_query == 'prot' and seq_type_db == 'DNA':
        blast_type = 'tblastn'
        call_blast(args, cat, 'tblastn')
        print ('Running tblastn...')
    else:
        print ("Error with sequence type, can't run blastn, blastp or tblastx")

    return blast_type

def call_blast(args, db, blast_type = 'blastn'):

    #cant use qcov_hsp_perc with cline (HSP = High- scoring Segment Pair )
    cmd = [blast_type,
    '-query', args.query,
    '-db', db,
    '-out', 'blast_out.txt',
    '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp',
    '-max_target_seqs', '500000',
    '-qcov_hsp_perc', args.percent_length,
    '-num_threads', args.threads]

    if blast_type == 'blastn':
        cmd += ['-perc_identity', args.percent_identity]
    call(cmd)

def tree(args, name, aln_suffix, rax = 'raxmlHPC-PTHREADS-SSE3', boots = '100', model = 'GTRGAMMA'):

    '''
    Makes a Raxml bipartitions tree
    Remove identical seqs from aln before running raxml
    Add them back to the final tree post raxml
    '''
    boots = str(boots)

    #set model. need to make this accessible one day, not hard coded...
    #Best tree
    if not os.path.exists(name + aln_suffix + '.reduced'):
        cmd = ['nice',
            rax,
            '-s', name + aln_suffix,
            '-n', name + '.tre',
            '-m', model,
            '-p', '12345',
            '-#', '20',
            '-T', args.threads]
        call(cmd)
    
    cmd = ['nice',
            rax,
            '-s', name + aln_suffix + '.reduced',
            '-n', name + '.tre',
            '-m', model,
            '-p', '12345',
            '-#', '20',
            '-T', args.threads]
    
    use_reduced = False
    if os.path.exists(name + aln_suffix + '.reduced'):
        if not os.path.exists(name + '_dup_info.txt'):#keep first one
            shutil.move('RAxML_info.' + name + '.tre', name + '_dup_info.txt')
        for rax_file in glob('RAxML_*'):
            os.remove(rax_file)
        call(cmd) #use the reduced aln so no dups
        print ('USING REDUCED!!!!!!!!!!!!!!!!')
        use_reduced = True
    
    if int(boots) > 0:

        if use_reduced:
            aln_4_boot = name + aln_suffix + '.reduced'
        else:
            aln_4_boot = name + aln_suffix
        
        #boots
        cmd = ['nice',
                rax,
                '-s', aln_4_boot,
                '-n', name + '_b.tre',
                '-m', model,
                '-p', '12345',
                '-b', '12345',
                '-#', boots,
                '-T', args.threads]
        call(cmd)
        #draw bipartitions on the best ML tree
        cmd = ['nice',
                rax,
                '-m', model,
                '-p', '12345',
                '-f', 'b',
                '-t', 'RAxML_bestTree.' + name + '.tre',
                '-z', 'RAxML_bootstrap.' + name + '_b.tre',
                '-n', name + '_bi.tre',
                '-T', args.threads]
        call(cmd)
        d = defaultdict(list)
        if os.path.exists(name + '_dup_info.txt'):
            with open(name + '_dup_info.txt','r')  as fin:
                for line in fin:
                    if line.startswith('IMPORTANT WARNING: Sequences'):
                        in_tree, _, not_int_tree = line.strip().split()[3:6]
                        d[in_tree].append(not_int_tree)
            for tre_file in ['RAxML_bestTree.' + name + '.tre',
                            'RAxML_bipartitionsBranchLabels.' + name + '_bi.tre',
                            'RAxML_bipartitions.' + name + '_bi.tre']:
                try:
                    with open(tre_file, 'r') as fin:
                        for line in fin:
                            tre = line.strip()
                            to_join = []
                            for sample in d:
                                to_join.append('(' + sample + ':0.0')
                                for sample2 in d.get(sample):
                                    to_join.append(',')
                                    to_join.append(sample2 + ':0.0')
                                    new = ''.join(to_join)+')'
                                    tre = tre.replace(sample, new)
                                    to_join = []
                    if tre_file == 'RAxML_bestTree.' + name + '.tre':
                        with open(name + '_' + tre_file.replace(name + '.tre', 'tre'), 'w') as fout:
                            fout.write(tre)
                    else:
                        with open(name + '_' + tre_file.replace(name + '_bi.tre', 'tre'), 'w') as fout:
                            fout.write(tre)
                except:
                    pass
        else:
            assert not use_reduced
            for tre_file in ['RAxML_bestTree.' + name + '.tre',
                            'RAxML_bipartitionsBranchLabels.' + name + '_bi.tre',
                            'RAxML_bipartitions.' + name + '_bi.tre']:
                if tre_file == 'RAxML_bestTree.' + name + '.tre':#copy instaed of move as that the output is same
                    shutil.copy(tre_file, name + '_' + tre_file.replace(name + '.tre', 'tre'))
                else:
                    shutil.copy(tre_file, name + '_' + tre_file.replace(name + '_bi.tre', 'tre'))
        #clean up
        if not os.path.exists('logs'):
            os.makedirs('logs')
        for log in glob('*.RUN.*'):
            shutil.move(log, 'logs/'+log)




if __name__ == "__main__":
        main()
