import csv


class Graph:

	'''
	idea for Graph class from http://code.activestate.com/recipes/569151-simple-graph-for-python/
	'''
	def __init__(self):
		self.node = {}
		self.neighbor = {}

		
	def add_node(self,node):
		self.node[node] = True


	def add_edge(self,node,nodesecond):
		try:
			self.neighbor[node].append(nodesecond)
		except:
			self.neighbor[node] = []
			self.neighbor[node].append(nodesecond)
		try:
			self.neighbor[nodesecond].append(node)
		except:
			self.neighbor[nodesecond] = []
			self.neighbor[nodesecond].append(node)


	def neighbor(self,node):
		return neighbor[node]

class Read_CSV:


	def __init__(self, inp_f, delimit):
		self.inp_f=inp_f
		self.delimit=delimit


	def read_csv(self):
		inp = open(self.inp_f, 'r')
		table = []
		for row in csv.reader(inp,delimiter=self.delimit):
			table.append(row)
		inp.close()
		return table

class Write_CSV:


	def __init__(self,table,output_file,delimiter):
		self.table=table
		self.output_file=output_file
		self.delimiter=delimiter


	def write_csv(self):
		output_file = open(self.output_file, 'w')
		writer = csv.writer(output_file,delimiter=self.delimiter)
		for row in self.table:
			writer.writerow(row)
		output_file.close()