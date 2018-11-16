from ROOT import *
from math import *
from array import array
import itertools
import numpy
import time


start_time = time.time()

f = TFile.Open("../data/TrainLinearData.root")
t = f.Get("DataTTree")
nentries = t.GetEntries()
windowsize = [1, 2, 3, 4]
addition_sp = 7
position = -1


class Event:
	def __init__(self, label, var1, var2):
		self.label = label
		self.var1  = var1
		self.var2  = var2

class BinPoint:
	def __init__(self, index, label, dist):
		self.index = index
		self.label = label
		self.dist  = dist


class NodeID:
	def __init__(self, var_type, index, bin_combination):
		self.varid = var_type
		self.index = index
		self.bins  = bin_combination


class ID3TreeNode:
	def __init__(self, label, var_type, index, posi, plane_para, ancestor_posi, offspring_list):
		self.label = label
		self.varid = var_type
		self.index = index # the id of the split bins = BinPoint.index = NodeID.index
		self.posi  = posi
		self.plane = plane_para
		self.ances = ancestor_posi
		self.offsp = offspring_list


def read_data():
	data_list = []
	for i in range(nentries):
		t.GetEntry(i)
		ivar1 = []
		for j in range(len(t.var1)):
			ivar1.append(t.var1[j])
		ivar2 = []
		for j in range(len(t.var2)):
			ivar2.append(t.var2[j])
		ilabel = t.label
		ievent = Event(ilabel, ivar1, ivar2)
		data_list.append(ievent)
	return data_list


def bin_combinations(var):
	n = len(var)
	bin_combinations_list = []
	for i in range(n):
		for j in windowsize:
			if i+j > n:
				break
			temp_list = []
			for k in range(j):
				temp_list.append(i+k)
			bin_combinations_list.append(temp_list)
	return bin_combinations_list


def find_centroid(events, label, current_events_list):
	var1_nbin = len(events[0].var1)
	var2_nbin = len(events[0].var2)
	centroid = Event(label, [0]*var1_nbin, [0]*var2_nbin)
	count = 0
	for i in current_events_list:
		if events[i].label != label:
			continue
		for j in range(var1_nbin):
			centroid.var1[j] += events[i].var1[j]
		for j in range(var2_nbin):
			centroid.var2[j] += events[i].var2[j]
		count += 1

	for j in range(var1_nbin):
		centroid.var1[j] /= count

	for j in range(var2_nbin):
		centroid.var2[j] /= count

	return centroid


#distance projection to the space given by bin_combination
def point_point_dist(var_a, var_b, bin_combination):
	distance = 0
	for i in bin_combination:
		distance += (var_a[i] - var_b[i])*(var_a[i] - var_b[i])
	return sqrt(distance)


#plane equation: [0]*x0 + [1]*x1 + ... + [n]*xn = 1
#this distance has sign, "+" means above the plane; "-" means below the plane
def point_plane_dist(var, plane_para):
	distance = 0
	for i in range(len(var)):
		distance += plane_para[i]*var[i]
	#distance = fabs(distance - 1)
	return distance-1


def sort_by_dist(bin_point):
	return fabs(bin_point.dist)


def find_split_points(events, bin_combination, current_events_list, var_type):
	centroids = [find_centroid(events, 0, current_events_list), find_centroid(events, 1, current_events_list)]
	bin_points = []
	for i in current_events_list:
		label = events[i].label
		centroid = centroids[1-label]
		if var_type == 1:
			var = events[i].var1
			cen = centroid.var1
		elif var_type == 2:
			var = events[i].var2
			cen = centroid.var2
		else:
			print("find_split_points(): no such variable!"),
			print(type(var_type))
			return None

		distance = point_point_dist(var, cen, bin_combination)
		bin_point = BinPoint(i, label, distance)
		bin_points.append(bin_point)

	bin_points.sort(key = sort_by_dist)
	sp = min(addition_sp + len(bin_combination), len(bin_points))
	split_points_list = []
	for i in range(sp):
		split_points_list.append(bin_points[i].index)

	return split_points_list
	

def cal_entropy(nsub, nall=1):
	if(nall == 0):
		return 0

	p = float(nsub)/float(nall)
	if(p > 1 or p < 0):
		print("cal_entropy(): wrong input probability!")
		return None
	if(p == 0 or p == 1):
		return 0
	else:
		H_entropy = - p*log(p) - (1-p)*log(1-p)
		return H_entropy/log(2)


def average_sub_entropy(events, para, current_events_list, var_type):
	current_size = len(current_events_list)
	sub_size = [0, 0]
	sub_label1_size = [0, 0]
	for i in current_events_list:
		if var_type == 1:
			var = events[i].var1
		elif var_type == 2:
			var = events[i].var2
		else:
			print("cal_entropy(): no such variable!")
			return None
		para_size = len(para)
		if para_size != len(var):
			print("cal_entropy(): para size is not consistent with varable size!"),
			print(len(var)),
			print(para_size)
			return None

		distance = point_plane_dist(var, para)
		label = events[i].label
		if distance < 0:
			sub_size[0] += 1
			if label == 1:
				sub_label1_size[0] += 1
		else:
			sub_size[1] += 1
			if label == 1:
				sub_label1_size[1] += 1

	H_average_sub_entropy = 1.0*sub_size[0]/current_size*cal_entropy(sub_label1_size[0], sub_size[0]) + sub_size[1]/current_size*cal_entropy(sub_label1_size[1], sub_size[1])

	return H_average_sub_entropy


def cal_split_plane(events, compact_split_points_list, bin_combination, var_type):
	size = len(compact_split_points_list)
	if size != len(bin_combination):
		print("cal_split_plane(): too much or too few points to calculate the plane!")
		return None
	points_matrix = []
	for i in compact_split_points_list:
		if var_type == 1:
			var = events[i].var1
		elif var_type == 2:
			var = events[i].var2
		else:
			print("cal_split_plane(): no such variable!")
			return None
		one_point = []
		for j in bin_combination:
			one_point.append(var[j])
		points_matrix.append(one_point)

	if size == 0:
		return [0]*len(var)

	if size == 1:
		plane_compact_para = [1.0/points_matrix[0][0]]
	else:
		try:
			inversed_points_matrix = numpy.linalg.inv(points_matrix)
		except:
			print("cal_split_plane(): matrix inversion problem!")
			return None
		RHS_ones = [1]*size
		plane_compact_para = numpy.matmul(inversed_points_matrix, RHS_ones)

	plane_para = [0]*len(var)
	for i in range(size):
		plane_para[bin_combination[i]] = plane_compact_para[i]

	return plane_para


def find_best_plane(events, split_points_list, bin_combination, current_events_list, var_type, entropy_cut = 1):
	best_entropy = min(1, entropy_cut) #the smaller the better
	best_plane_para = [0]*(var_type+3)
	size = len(bin_combination)
	for compact_split_points_list in itertools.combinations(split_points_list, size):
		current_plane_para = cal_split_plane(events, compact_split_points_list, bin_combination, var_type)
		current_entropy    = average_sub_entropy(events, current_plane_para, current_events_list, var_type)
		if current_entropy <= best_entropy:
			best_entropy    = current_entropy
			best_plane_para = current_plane_para

	#print(best_entropy)
	return best_plane_para


def refine_split_points(events, bin_combination, current_events_list, var_type, plane_para):
	bin_points = []
	for i in current_events_list:
		label = events[i].label
		if var_type == 1:
			var = events[i].var1
		elif var_type == 2:
			var = events[i].var2
		else:
			print("find_split_points(): no such variable!"),
			print(type(var_type))
			return None

		distance = abs(point_plane_dist(var, plane_para))
		bin_point = BinPoint(i, label, distance)
		bin_points.append(bin_point)

	bin_points.sort(key = sort_by_dist)
	sp = min(addition_sp + len(bin_combination), len(bin_points))
	split_points_list = []
	for i in range(sp):
		split_points_list.append(bin_points[i].index)

	return split_points_list


def refine_best_plane(events, bin_combination, current_events_list, var_type, best_plane_para):
	best_entropy = average_sub_entropy(events, best_plane_para, current_events_list, var_type)
	best_split_points = refine_split_points(events, bin_combination, current_events_list, var_type, best_plane_para)
	return find_best_plane(events, best_split_points, bin_combination, current_events_list, var_type, best_entropy)


def find_nodes_list(events):
	nodes_list = []
	node_id = 0
	var1_bin_combinations_list = bin_combinations(events[0].var1)
	var2_bin_combinations_list = bin_combinations(events[0].var2)
	var_type = 1
	for i in var1_bin_combinations_list:
		one_node = NodeID(var_type, node_id, i)
		nodes_list.append(one_node)
		node_id += 1
	var_type = 2
	for i in var2_bin_combinations_list:
		one_node = NodeID(var_type, node_id, i)
		nodes_list.append(one_node)
		node_id += 1	

	return nodes_list


def find_best_split_node(events, ances_node_posi, current_events_list, current_nodes_list):
	best_entropy = 1 #the smaller the better
	best_plane   = []
	best_nodeid  = -1
	best_varid   = 0
	for trial_node in current_nodes_list:
		trial_split_points = find_split_points(events, trial_node.bins, current_events_list, trial_node.varid) 
		trial_best_plane = find_best_plane(events, trial_split_points, trial_node.bins, current_events_list, trial_node.varid)
		for n_refine in range(3):
			trial_best_plane = refine_best_plane(events, trial_node.bins, current_events_list, trial_node.varid, trial_best_plane)
 
		trial_entropy    = average_sub_entropy(events, trial_best_plane, current_events_list, trial_node.varid)
		if(trial_entropy < best_entropy):
			best_plane   = trial_best_plane
			best_entropy = trial_entropy
			best_nodeid  = trial_node.index
			best_varid   = trial_node.varid

	global position
	best_label = find_majority_label(events, current_events_list)
	best_split_node = ID3TreeNode(best_label, best_varid, best_nodeid, position, best_plane, ances_node_posi, [])

	return best_split_node	


def find_majority_label(events, current_events_list):
	label0 = 0
	label1 = 0
	for i in current_events_list:
		if events[i].label == 0:
			label0 += 1
		if events[i].label == 1:
			label1 += 1
	if label0 > label1:
		return 0
	else:
		return 1


def list_is_pure(events, current_events_list):
	n_label0 = 0
	n_label1 = 0
	for i in current_events_list:
		if events[i].label == 0:
			n_label0 += 1
		if events[i].label == 1:
			n_label1 += 1
	if n_label0 == 0 or n_label1 == 0:
		return True
	else:
		return False



def build_ID3_tree(events, nodes_storage, ances_node, current_events_list, current_nodes_list):
	global position

	if nodes_storage == []:
		position += 1
		ances_node = find_best_split_node(events, -1, current_events_list, current_nodes_list)
		nodes_storage.append(ances_node)
	
	if len(current_events_list) <= 5 or len(current_nodes_list) == 0 or list_is_pure(events, current_events_list):
		position += 1
		terminal_label = find_majority_label(events, current_events_list)
		terminal_node  = ID3TreeNode(terminal_label, ances_node.varid, -2, position, [], ances_node.posi, []) #index = -2 means the termial node
		nodes_storage.append(terminal_node)
		return None

	left_events_list = []
	right_events_list= []

	for i in current_events_list:
		if ances_node.varid == 1:
			var = events[i].var1
		elif ances_node.varid == 2:
			var = events[i].var2
		else:
			print("build_ID3_tree(): no such variable!")
			return None

		distance = point_plane_dist(var, ances_node.plane)
		if distance < 0:
			left_events_list.append(i)
		else:
			right_events_list.append(i)

	if list_is_pure(events, left_events_list):
		position += 1
		terminal_label = find_majority_label(events, left_events_list)
		terminal_node  = ID3TreeNode(terminal_label, ances_node.varid, -2, position, [], ances_node.posi, [])
		nodes_storage.append(terminal_node)
	else:
		position += 1
		left_best_split_node  = find_best_split_node(events, ances_node.posi, left_events_list, current_nodes_list) 
		nodes_storage.append(left_best_split_node)
		left_nodes_list = list(current_nodes_list)
		#count = 0
		#for one_node in current_nodes_list:
		#	if one_node.index == left_best_split_node.index:
		#		del left_nodes_list[count]
		#		break
		#	count += 1
		build_ID3_tree(events, nodes_storage, left_best_split_node, left_events_list, left_nodes_list)

	if list_is_pure(events, right_events_list):
		position += 1
		terminal_label = find_majority_label(events, right_events_list)
		terminal_node  = ID3TreeNode(terminal_label, ances_node.varid, -2, position, [], ances_node.posi, [])
		nodes_storage.append(terminal_node)
	else:
		position += 1
		right_best_split_node = find_best_split_node(events, ances_node.posi, right_events_list,current_nodes_list)
		nodes_storage.append(right_best_split_node)
		right_nodes_list = list(current_nodes_list)
		#count = 0
		#for one_node in current_nodes_list:
		#	if one_node.index == right_best_split_node.index:
		#		del right_nodes_list[count]
		#		break
		#	count += 1
		build_ID3_tree(events, nodes_storage, right_best_split_node, right_events_list, right_nodes_list)

	#print(position)
	return None


def main():
	events = read_data()
	full_events_list = range(nentries)
	full_nodes_list = find_nodes_list(events)
	tree_nodes = []
	root_node_index = -1
	build_ID3_tree(events, tree_nodes, root_node_index, full_events_list, full_nodes_list)

	for one_node in tree_nodes:
		#print(one_node.posi)
		tree_nodes[one_node.ances].offsp.append(one_node.posi)

	#save the ID3 tree to a TTree
	node_label    = array('i', [0])
	node_var_type = array('i', [0])
	node_index    = array('i', [0])
	node_posi     = array('i', [0])
	node_ances    = array('i', [0])
	node_n_offsp  = array('i', [0])
	node_offsp    = array('i', [0]*2)
	node_plane    = array('d', [0.0]*5)
	f_tree = TFile("ID3_tree.root", "recreate")
	t_tree = TTree("ID3TTree", "the pre-pruned the ID3 tree")
	t_tree.Branch("label"   , node_label   , "label/I")
	t_tree.Branch("var_type", node_var_type, "var_type/I")
	t_tree.Branch("index"   , node_index   , "index/I")
	t_tree.Branch("posi"    , node_posi    , "posi/I")
	t_tree.Branch("ances"   , node_ances   , "ances_position/I")
	t_tree.Branch("n_offsp" , node_n_offsp , "n_offsp/I")
	t_tree.Branch("offsp"   , node_offsp   , "offsprings[n_offsp]/I")
	t_tree.Branch("plane"   , node_plane   , "plane[5]/D")

	for one_node in tree_nodes:
		node_label[0]    = one_node.label
		node_var_type[0] = one_node.varid
		node_index[0]    = one_node.index
		node_posi[0]     = one_node.posi
		node_ances[0]    = one_node.ances
		node_n_offsp[0]  = len(one_node.offsp)
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


main()

print("--- %s seconds ---" % (time.time() - start_time))




