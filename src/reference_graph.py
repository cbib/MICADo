#!/usr/bin/env python
# coding=utf8
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import networkx as nx

####################################################
####################################################
## Parameters loading 
# Data directories
data_directory = "data/reference"