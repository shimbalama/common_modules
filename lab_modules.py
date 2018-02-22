#!/usr/bin/env python3

def main ():
        pass
def tree(args, model = 'GTRGAMMA'):

        '''
        Makes a Raxml bipartitions tree
	Remove identical seqs from aln before running raxml
	Add them back to the final tree post raxml
        '''
        #set model. need to make this accessible one day, not hard coded...
	#Best tree
        cmd = ['nice',
                'raxmlHPC-PTHREADS-AVX2',
                '-s', args.out+'_trimmed.aln',
                '-n', args.out+'.tre',
                '-m', model,
                '-p', '12345',
                '-#', '20',
                '-T', args.threads]
        call(cmd)
        cmd = ['nice',
                'raxmlHPC-PTHREADS-AVX2',
                '-s', args.out+'_trimmed.aln.reduced',
                '-n', args.out+'.tre',
                '-m', model,
                '-p', '12345',
                '-#', '20',
                '-T', args.threads]
        if os.path.exists(args.out + '_trimmed.aln.reduced'):
                shutil.move('RAxML_info.' + args.out + '.tre', 'dup_info.txt')
                for rax in glob('RAxML_*'):
                        os.remove(rax)
                call(cmd) #use the reduced aln so no dups
                print ('USING REDUCED!!!!!!!!!!!!!!!!')

        #boots
        cmd = ['nice',
                'raxmlHPC-PTHREADS-AVX2',
                '-s', args.out+'_trimmed.aln',
                '-n', args.out+'_boot.tre',
                '-m', model,
                '-p', '12345',
                '-b', '12345',
                '-#', '100',
                '-T', args.threads]
        call(cmd)
        #draw bipartitions on the best ML tree
        cmd = ['nice',
                'raxmlHPC-PTHREADS-AVX2',
                '-m', model,
                '-p', '12345',
                '-f', 'b',
                '-t', 'RAxML_bestTree.' + args.out + '.tre',
                '-z', 'RAxML_bootstrap.' + args.out + '_boot.tre',
                '-n', args.out+'_boots_on.tre',
                '-T', args.threads]
        call(cmd)

        d=collections.defaultdict(list)
        with open('dup_info.txt','r')  as fin:
                for line in fin:
                        if line.startswith('IMPORTANT WARNING: Sequences'):
                                in_tree, _, not_int_tree = line.strip().split()[3:6]
                                d[in_tree].append(not_int_tree)
        for tre_file in ['RAxML_bestTree.'+args.out+'.tre',
                        'RAxML_bestTree_bipartitionsBranchLabels.'+args.out+'_boots_on.tre',
                        'RAxML_bestTree_bipartitions.'+args.out+'_boots_on.tre']:
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
                        with open(tre_file + '.all.tre', 'w') as fout:
                                fout.write(tre)
                except:
                        pass

if __name__ == "__main__":
        main()
