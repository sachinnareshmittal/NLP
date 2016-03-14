from collections import defaultdict
import sys
import numpy

class HMM():
	def __init__(self,train,test):
		self.ftrain = train
		self.ftest = test


	def count(self):
		self.bigram = defaultdict(int)
		self.unigram = defaultdict(int)
		self.trigram = defaultdict(int)

		p_previous_tag = START
		previous_tag = START
		current_tag = START


		for line in open(self.ftrain,'r'):
			line = line.strip()
			p_previous_tag = START
			previous_tag = START
			current_tag = START

			for l in line.split(" "):
				l = l.strip()
				word = ''
				tag = ''


