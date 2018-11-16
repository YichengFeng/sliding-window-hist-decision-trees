from ROOT import *
from math import *
from array import array
import itertools
import numpy


f = TFile.Open("../../../data/PruneNonlinearData.root")
t = f.Get("DataTTree")
nentries = t.GetEntries()

f_ID3 = TFile.Open("../train/ID3_tree.root")
t_ID3 = f_ID3.Get("ID3TTree")
n_nodes = t_ID3.GetEntries()
var1_nbin = 4
var2_nbin = 5


class Event:
	def __init__(self, label, var1, var2):
		self.label = label
		self.var1  = var1
		self.var2  = var2

class NodeID:
	def __init__(self, var_type, index, bin_combination):
		self.varid = var_type
		self.index = index
		self.bins  = bin_combination


class ID3TreeNodePrune:
	def __init__(self, label, var_type, index, posi, plane_para, ancestor_posi, offspring_list, N_fraction, accuracy, best_offsp_weighted_accu):
		self.label = label
		self.varid = var_type
		self.index = index # the id of the split bins = BinPoint.index = NodeID.index
		self.posi  = posi
		self.plane = plane_para
		self.ances = ancestor_posi
		self.offsp = offspring_list
		self.Nfrac = N_fraction
		self.accu  = accuracy
		self.best_offsp_weighted_accu = best_offsp_weighted_accu


def read_data():
	data_list = []
	for i in range(nentries):
		t.GetEntry(i)
		ivar1 = []
		for j in range(len(t.var1)):
			ivar1.append(t.var1[j])
		ivar2 = []
		#for j in range(len(t.var2)):
		#	ivar2.append(t.var2[j])
		ilabel = t.label
		ievent = Event(ilabel, ivar1, ivar2)
		data_list.append(ievent)
	return data_list


def read_ID3_tree():
	ID3_tree_prune = []
	for i in range(n_nodes):
		t_ID3.GetEntry(i)
		label    = t_ID3.label
		var_type = t_ID3.var_type
		index    = t_ID3.index
		posi     = t_ID3.posi
		ances    = t_ID3.ances
		n_offsp  = t_ID3.n_offsp
		offsp    = []
		plane    = []
		if var_type == 1:
			var_nbin = var1_nbin
		elif var_type == 2:
			var_nbin = var2_nbin
		else:
			print("read_ID3_tree(): no such variable!")
			print(var_type)
			return None
		for j in range(n_offsp):
			offsp.append(t_ID3.offsp[j])
		for j in range(var_nbin):
			plane.append(t_ID3.plane[j])
		ID3_tree_prune.append(ID3TreeNodePrune(label, var_type, index, posi, plane, ances, offsp, 0, 0, 0))

	return ID3_tree_prune


def point_plane_dist(var, plane_para):
	distance = 0
	for i in range(len(var)):
		distance += plane_para[i]*var[i]
	#distance = fabs(distance - 1)
	return distance-1


def cal_frac_accu(events, current_events_list, old_tree, node_position):
	current_amount = len(current_events_list)
	current_node = old_tree[node_position]
	node_label = current_node.label
	current_node.Nfrac = current_amount
	if current_amount == 0:
		accuracy = 0 #float("nan")
	else:
		correct_prediction = 0
		for i in current_events_list:
			if events[i].label == node_label:
				correct_prediction += 1
		accuracy = 1.0*correct_prediction/current_amount
	current_node.accu = accuracy
	current_node.best_offsp_weighted_accu = 0 

	if current_node.index == -2:
		return None

	left_events_list = []
	right_events_list= []

	for i in current_events_list:
		var_type = current_node.varid
		if var_type == 1:
			var = events[i].var1
		elif var_type == 2:
			var = events[i].var2
		else:
			print("cal_frac_accu(): no such variable!")
			return None

		distance = point_plane_dist(var, current_node.plane)
		if distance < 0:
			left_events_list.append(i)
		else:
			right_events_list.append(i)

	if len(current_node.offsp) == 0:
		print("cal_frac_accu(): no offspring!")
		return None
	if len(current_node.offsp) >= 1:
		cal_frac_accu(events, left_events_list, old_tree, current_node.offsp[0])
	if len(current_node.offsp) >= 2:
		cal_frac_accu(events, right_events_list,old_tree, current_node.offsp[1])

	return None


def find_best_offsp_weighted_accu(old_tree, node_position, prune_node_list):
	current_node = old_tree[node_position]
	weighted_accu = current_node.accu*current_node.Nfrac
	if current_node.index == -2:
		current_node.best_offsp_weighted_accu = weighted_accu
		return weighted_accu

	left_weighted_accu = 0
	right_weighted_accu= 0

	if len(current_node.offsp) == 0:
		print("find_best_offsp_weighted_accu(): no offspring!")
		return weighted_accu
	if len(current_node.offsp) >= 1:
		left_weighted_accu = find_best_offsp_weighted_accu(old_tree, current_node.offsp[0], prune_node_list)
	if len(current_node.offsp) >= 2:
		right_weighted_accu= find_best_offsp_weighted_accu(old_tree, current_node.offsp[1], prune_node_list)

	#weighted_accu = max(weighted_accu, (left_weighted_accu+right_weighted_accu))
	if (left_weighted_accu + right_weighted_accu) > weighted_accu:
		weighted_accu = left_weighted_accu + right_weighted_accu
	else:
		prune_node_list.append(node_position)

	current_node.best_offsp_weighted_accu = weighted_accu

	return weighted_accu


def prune_sub_tree(old_tree, node_position):
	current_node = old_tree[node_position]
	current_node.index = -2
	for offsp_posi in current_node.offsp:
		if offsp_posi <1:
			continue
		prune_sub_tree(old_tree, offsp_posi)

	return None


def find_prune_nodes(old_tree, node_position):
	current_node = old_tree[node_position]
	if(current_node.index == -2):
		return None

	weighted_accu = current_node.accu*current_node.Nfrac
	if weighted_accu >= current_node.best_offsp_weighted_accu:
		prune_sub_tree(old_tree, node_position)

	return None
		

def prune_main():
	events = read_data()
	old_tree = read_ID3_tree()
	fulllist = range(nentries)

	prune_posi_list = []
	cal_frac_accu(events, fulllist, old_tree, 0)
	find_best_offsp_weighted_accu(old_tree, 0, prune_posi_list)
	for prune_posi in prune_posi_list:
		prune_sub_tree(old_tree, prune_posi)

	#find_prune_nodes(old_tree, 0)

	#save the ID3 tree to a TTree
	node_label    = array('i', [0])
	node_var_type = array('i', [0])
	node_index    = array('i', [0])
	node_posi     = array('i', [0])
	node_ances    = array('i', [0])
	node_n_offsp  = array('i', [0])
	node_n_events = array('i', [0])
	node_accu     = array('d', [0.0])
	node_best_offsp_weighted_accu = array('d', [0.0])
	node_offsp    = array('i', [0]*2)
	node_plane    = array('d', [0.0]*5)
	f_tree = TFile("postpruned_ID3_tree.root", "recreate")
	t_tree = TTree("PrunedID3TTree", "the pre-pruned the ID3 tree")
	t_tree.Branch("label"   , node_label   , "label/I")
	t_tree.Branch("var_type", node_var_type, "var_type/I")
	t_tree.Branch("index"   , node_index   , "index/I")
	t_tree.Branch("posi"    , node_posi    , "posi/I")
	t_tree.Branch("ances"   , node_ances   , "ances_position/I")
	t_tree.Branch("n_offsp" , node_n_offsp , "n_offsp/I")
	t_tree.Branch("accu"    , node_accu    , "accuracy/D")
	t_tree.Branch("n_events", node_n_events, "n_events/I")
	t_tree.Branch("right_decision", node_best_offsp_weighted_accu, "right_decision/D")
	t_tree.Branch("offsp"   , node_offsp   , "offsprings[n_offsp]/I")
	t_tree.Branch("plane"   , node_plane   , "plane[5]/D")

	for one_node in old_tree:
		node_label[0]    = one_node.label
		node_var_type[0] = one_node.varid
		node_index[0]    = one_node.index
		node_posi[0]     = one_node.posi
		node_ances[0]    = one_node.ances
		node_n_offsp[0]  = len(one_node.offsp)
		node_accu[0]     = one_node.accu
		node_n_events[0] = one_node.Nfrac
		node_best_offsp_weighted_accu[0] = one_node.best_offsp_weighted_accu
		for i in range(2):
			node_offsp[i] = 0
		for i in range(5):
			node_plane[i] = 0
		for i in range(node_n_offsp[0]):
			node_offsp[i] = one_node.offsp[i]
		for i in range(len(one_node.plane)):
			node_plane[i] = one_node.plane[i]
		t_tree.Fill()

	f_tree.Write()
	f_tree.Close()
	
	return

prune_main()

 


##############################
##test
#events = read_data()
#old_tree = read_ID3_tree()
#print(old_tree[0].plane)

#fulllist = range(nentries)
#cal_frac_accu(events, fulllist, old_tree, 0)
#print(old_tree[70].accu)
