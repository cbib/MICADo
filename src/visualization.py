#!/usr/bin/env python
# coding=utf8
import networkx as nx
import itertools

import logging

if "logger" not in globals():
	logger = logging.getLogger('GRAPHMERGER')
	logger.setLevel(logging.DEBUG)
	ch = logging.StreamHandler()
	ch.setLevel(logging.DEBUG)
	formatter = logging.Formatter('%(asctime)s - %(filename)s - %(message)s', "%Y-%m-%d %H:%M:%S")
	ch.setFormatter(formatter)
	logger.addHandler(ch)

logger = logging.getLogger('GRAPHMERGER')


# Merge reference graph
def merge_reference_graph(reference_graph):
	# create a copy 
	merged_graph = nx.DiGraph(reference_graph)
	# singletons ? 
	out_d = merged_graph.out_degree()
	in_d = merged_graph.in_degree()
	singletons = [x for x, deg in in_d.items() if (deg == 1) and out_d[x] == 1]
	# bunch them
	induced = merged_graph.subgraph(singletons).to_undirected()
	buckets = itertools.ifilter(lambda x: len(x) > 1, nx.networkx.connected_components(induced.to_undirected()))
	is_entry_point = lambda x: reference_graph.pred[x].items()[0][0] not in singletons
	is_exit_point = lambda x: reference_graph.succ[x].items()[0][0] not in singletons
	n_buckets = 0
	# Compact them in the merged graph
	for bunch in buckets:
		entry_point = reference_graph.predecessors_iter(
			itertools.ifilter(is_entry_point, bunch).next()).next()  # we know there's exactly 1 predecessor and thus 1 entry_point whose pred is in the graph
		exit_point = reference_graph.successors_iter(
			itertools.ifilter(is_exit_point, bunch).next()).next()  # we know there's exactly 1 successor and thus 1 exit point  whose succ exists in the graph
		# print "Compacting chain (unordered): %s"%" / ".join(bunch)
		meta_node_lbl = "_".join(bunch)
		# color the meta node with the 'ref_list' attribut from the first node of the meta node
		ref_list_meta_node = []
		fonction_ref = lambda x: merged_graph.node[x]['ref_list']
		for node in bunch:
			ref_list_meta_node += fonction_ref(node)
		merged_graph.add_node(meta_node_lbl, length=len(bunch), ref_list=set(ref_list_meta_node))  # je donne un poids au noeud
		merged_graph.remove_nodes_from(bunch)  # j'efface tout
		merged_graph.add_edge(entry_point, meta_node_lbl)  # j'ajoute une arrete entre mon noeud précédent le méta noeud et le méta noeud
		merged_graph.add_edge(meta_node_lbl, exit_point)  # j'ajoute une arrete entre mon noeuds suivant le méta noeud le méta noeud
		n_buckets += 1
	logger.info("Found and compacted %d linear chains", n_buckets)
	logger.info("Reducing graph size from %d to %d", len(reference_graph), len(merged_graph))
	return (merged_graph)


# Format reference graph for cytoscape visualization
def reference_graph_visualization_formatting(reference_graph):
	graph2visu = nx.DiGraph(reference_graph)
	for node in graph2visu.nodes():
		graph2visu.node[node]['ref_list'] = str(";".join(sorted(graph2visu.node[node]['ref_list'])))
	return (graph2visu)

# Format reference graph merged for cytoscape visualization
def reference_graph_merged_visualization_formatting(reference_graph):
	graph2visu = nx.DiGraph(reference_graph)
	for node in graph2visu.nodes():
		graph2visu.node[node]['ref_list'] = str(";".join(sorted(graph2visu.node[node]['ref_list'])))
		if graph2visu.node[node].get('length', 0) == 0:
			graph2visu.node[node]['length'] = 1
	return (graph2visu)

# Merge individu graph and color(attribut = ref_list) nodes
def merge_individu_graph(input_graph, reference_graph):
	# create a copy 
	merged_graph = nx.DiGraph(input_graph)
	# singletons ? 
	out_d = merged_graph.out_degree()
	in_d = merged_graph.in_degree()
	singletons = [x for x, deg in in_d.items() if (deg == 1) and out_d[x] == 1]
	# bunch them
	induced = merged_graph.subgraph(singletons).to_undirected()
	buckets = itertools.ifilter(lambda x: len(x) > 1, nx.networkx.connected_components(induced.to_undirected()))
	is_entry_point = lambda x: input_graph.pred[x].items()[0][0] not in singletons
	is_exit_point = lambda x: input_graph.succ[x].items()[0][0] not in singletons
	n_buckets = 0
	# Compact them in the merged graph
	for bunch in buckets:
		entry_point = input_graph.predecessors_iter(
			itertools.ifilter(is_entry_point, bunch).next()).next()  # we know there's exactly 1 predecessor and thus 1 entry_point whose pred is in the graph
		exit_point = input_graph.successors_iter(
			itertools.ifilter(is_exit_point, bunch).next()).next()  # we know there's exactly 1 successor and thus 1 exit point  whose succ exists in the graph
		# print "Compacting chain (unordered): %s"%" / ".join(bunch)
		meta_node_lbl = "_".join(bunch)
		# color the meta node with the 'ref_list' attribut from all nodes present of the meta node
		read_list_meta_node = []
		fonction_read_list = lambda x: merged_graph.node[x]['read_list_n']
		ref_list_meta_node = []
		fonction_ref = lambda x: reference_graph.node[x]['ref_list']
		condition = 0
		for node in bunch:
			read_list_meta_node += fonction_read_list(node)
			if node in reference_graph:
				condition = 1
				ref_list_meta_node += fonction_ref(node)
		if condition == 1:
			merged_graph.add_node(meta_node_lbl, length=len(bunch), ref_list=set(ref_list_meta_node), read_list_n=len(set(read_list_meta_node)))  # je donne un poids au noeud
		else:
			merged_graph.add_node(meta_node_lbl, length=len(bunch), ref_list=set("alt"), read_list_n=len(set(read_list_meta_node)))  # je donne un poids au noeud
		merged_graph.remove_nodes_from(bunch)  # j'efface tout
		merged_graph.add_edge(entry_point, meta_node_lbl)  # j'ajoute une arrete entre mon noeud précédent le méta noeud et le méta noeud
		merged_graph.add_edge(meta_node_lbl, exit_point)  # j'ajoute une arrete entre mon noeuds suivant le méta noeud le méta noeud
		n_buckets += 1
	logger.info("Found and compacted %d linear chains",n_buckets)
	logger.info("Reducing graph size from %d to %d",len(input_graph), len(merged_graph))
	return merged_graph

# Format individu graph for cytoscape visualization
def individu_graph_visualization_formating(input_graph, reference_graph):
	graph2visu = nx.DiGraph(input_graph)
	for node in graph2visu.nodes():
		if node in reference_graph:
			graph2visu.node[node]['ref_list'] = str(";".join(sorted(reference_graph.node[node]['ref_list'])))
		else:
			graph2visu.node[node]['ref_list'] = "alt"
		graph2visu.node[node]['read_list_n'] = len(set(graph2visu.node[node]['read_list_n']))
		del graph2visu.node[node]['fastq_id']
	return graph2visu

# Format individu graph for cytoscape visualization
def individu_graph_merged_visualization_formating(input_graph, reference_graph):
	graph2visu = nx.DiGraph(input_graph)
	for node in graph2visu.nodes():
		if graph2visu.node[node].get('ref_list', 0) != 0:
			graph2visu.node[node]['ref_list'] = str(";".join(sorted(graph2visu.node[node]['ref_list'])))
		else:
			graph2visu.node[node]['read_list_n'] = len(set(graph2visu.node[node]['read_list_n']))
			graph2visu.node[node]['length'] = 1
			del graph2visu.node[node]['fastq_id']
			if node in reference_graph:
				graph2visu.node[node]['ref_list'] = str(";".join(sorted(reference_graph.node[node]['ref_list'])))
			else:
				graph2visu.node[node]['ref_list'] = "alt"
	return graph2visu
