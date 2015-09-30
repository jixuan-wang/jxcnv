import argparse
import numpy as np

def svd(args):
	filename = args.datafile 
	
	# count the columns number of the data file
	f = open(filename)
	line = f.readline()
	colsnum = len(line.split('\t'))

	# skip 1st row and 1st column
	data = np.loadtxt(filename, dtype=np.float, delimiter='\t', skiprows=1, usecols=range(1, colsnum)) 
	# test for conifer
	# data = np.loadtxt(filename, dtype=np.float, delimiter='\t')

	### for test
	data = np.transpose(data) 

	# svd transform
	# comp_removed = args.svd
	U, S, Vt = np.linalg.svd(data, full_matrices=False)
	
	# new_S = np.diag(np.hstack([np.zeros([comp_removed]), S[comp_removed:]]))
	index = S < 0.7 * np.mean(S)
	new_S = np.diag(S * index)
	
	# reconstruct data matrix
	data = np.dot(U, np.dot(new_S, Vt))
	   
	# save to files
	file_u = open(args.output + '.U', 'w')
	file_s = open(args.output + '.S', 'w')
	file_vt = open(args.output + '.Vt', 'w')
	file_svd = open(args.output + '.SVD', 'w')

	np.savetxt(file_u, U, delimiter='\t')
	np.savetxt(file_s, S, delimiter='\t')
	np.savetxt(file_vt, Vt, delimiter='\t')
	np.savetxt(file_svd, np.transpose(data), delimiter='\t')


	file_u.close()
	file_s.close()
	file_vt.close()
	file_svd.close()	
	

parser = argparse.ArgumentParser(prog='jxcnv', description='Designed by jx.')
subparsers = parser.add_subparsers()

#SVD
svd_parser = subparsers.add_parser('svd', help="SVD")
svd_parser.add_argument('--datafile', required=True, help='')
svd_parser.add_argument('--output', required=True, help='')
# svd_parser.add_argument('--svd', type=int, required=True, help='Number of components to remove')
svd_parser.set_defaults(func=svd)

args = parser.parse_args()
args.func(args)