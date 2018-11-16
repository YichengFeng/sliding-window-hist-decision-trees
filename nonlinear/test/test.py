from ROOT import *
from math import *
from array import array
import itertools
import numpy


f = TFile.Open("../data/TestNonlinearData.root")
t = f.Get("DataTTree")
nentries = t.GetEntries()

f_ID3 = TFile.Open("../prune/postpruned_ID3_tree.root")
t_ID3 = f_ID3.Get("PrunedID3TTree")
#f_ID3 = TFile.Open("../train/ID3_tree.root")
#t_ID3 = f_ID3.Get("ID3TTree")
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
		#the last 3 features are not used here.


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


#the following function is not used here:
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
	current_node.best_offsp_weighted_accu = accuracy*current_amount

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


def find_best_offsp_weighted_accu(old_tree, node_position):
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
		left_weighted_accu = find_best_offsp_weighted_accu(old_tree, current_node.offsp[0])
	if len(current_node.offsp) >= 2:
		right_weighted_accu= find_best_offsp_weighted_accu(old_tree, current_node.offsp[1])

	#weighted_accu = max(weighted_accu, (left_weighted_accu+right_weighted_accu))
	if (left_weighted_accu + right_weighted_accu) > weighted_accu:
		weighted_accu = left_weighted_accu + right_weighted_accu
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
	if current_node.index == -2:
		return None

	weighted_accu = current_node.accu*current_node.Nfrac
	if weighted_accu >= current_node.best_offsp_weighted_accu:
		prune_sub_tree(old_tree, node_position)

	return None
		

def find_event_label(event, ID3_tree, node_position):
	current_node = ID3_tree[node_position]
	if current_node.index == -2:
		return current_node.label

	var_type = current_node.varid
	if var_type == 1:
		var = event.var1
	elif var_type == 2:
		var = event.var2
	else:
		print("cal_frac_accu(): no such variable!")
		return None

	if len(current_node.offsp) == 1:
		return find_event_label(event, ID3_tree, current_node.offsp[0])

	distance = point_plane_dist(var, current_node.plane)
	if distance < 0:
		return find_event_label(event, ID3_tree, current_node.offsp[0])
	else:
		return find_event_label(event, ID3_tree, current_node.offsp[1])



def cal_divider(ID3_tree, canvas):
	dividers = []
	canvas.cd()
	for node in ID3_tree:
		if node.index == -2:
			continue
		if node.varid == 2:
			continue
		plane = node.plane

		skip = False
		for i in range(2, var1_nbin):
			if plane[i] != 0:
				skip = True
		if skip:
			continue

		if plane[0] == 0 and plane[1] == 0:
			continue
		if plane[0] != 0 and plane[1] == 0:
			x = 1.0/plane[0]
			line = TLine(x, 0, x, 1)
		if plane[0] == 0 and plane[1] != 0:
			y = 1.0/plane[1]
			line = TLine(0, y, 1, y)
		if plane[0] != 0 and plane[1] != 0:
			x1 = 1.0/plane[0]
			x2 = (1.0 - plane[1])/plane[0]
			line = TLine(x1, 0, x2, 1)

			line.SetLineWidth(1)
			line.Draw('l')
			dividers.append(line)

	print('projection lines: '),
	print(len(dividers))

	return dividers



def main():
	events = read_data()
	ID3_tree = read_ID3_tree()
	fulllist = range(nentries)

	n_split_node = 0
	for node in ID3_tree:
		if node.index != -2:
			n_split_node += 1
	print('split nodes of postpruned_ID3_tree: '),
	print(n_split_node)

	hvar1bin12label0 = TGraph();
	hvar1bin12label0.SetMarkerStyle(21)
	hvar1bin12label0.SetMarkerSize(0.5)
	hvar1bin12label0.SetMarkerColor(kBlue)
	hvar1bin12label1 = TGraph();
	hvar1bin12label1.SetMarkerStyle(21)
	hvar1bin12label1.SetMarkerSize(0.5)
	hvar1bin12label1.SetMarkerColor(kRed)
	hvar1bin12wrong  = TGraph();
	hvar1bin12wrong.SetMarkerStyle(21)
	hvar1bin12wrong.SetMarkerSize(0.5)
	hvar1bin12wrong.SetMarkerColor(kGreen)

	#cal_frac_accu(events, fulllist, old_tree, 0)
	#find_best_offsp_weighted_accu(old_tree, 0)
	#find_prune_nodes(old_tree, 0)

	mistake0as1 = 0
	mistake1as0 = 0

	for event in events:
		classifier = find_event_label(event, ID3_tree, 0)
		if classifier == event.label:
			if classifier == 0:
				hvar1bin12label0.SetPoint(hvar1bin12label0.GetN(), event.var1[0], event.var1[1])
			if classifier == 1:
				hvar1bin12label1.SetPoint(hvar1bin12label1.GetN(), event.var1[0], event.var1[1])
		else:
			if classifier == 0:
				mistake0as1 += 1
			if classifier == 1:
				mistake1as0 += 1
			hvar1bin12wrong.SetPoint(hvar1bin12wrong.GetN(), event.var1[0], event.var1[1])

	cvar1 = TCanvas("var1bin12", "var1bin12", 600, 600)
	hvar1bin12label0.GetXaxis().SetRangeUser(0,1)
	hvar1bin12label0.GetYaxis().SetRangeUser(0,1)
	hvar1bin12label0.GetXaxis().SetTitle("bin 1 of variable 1 (V1_{1})")
	hvar1bin12label0.GetYaxis().SetTitle("bin 2 of variable 1 (V1_{2})")
	hvar1bin12label0.Draw("AP")
        hvar1bin12label1.Draw("P same")
        hvar1bin12wrong.Draw("P same")

	leg = TLegend(0.55, 0.7, 0.95, 0.95)
	leg.SetFillStyle(0)
	leg.SetBorderSize(0)
	leg.AddEntry(hvar1bin12label0, "#splitline{correct classification}{ with label 0}", 'p')
	leg.AddEntry(hvar1bin12label1, "#splitline{correct classification}{ with label 1}", 'p')
	leg.AddEntry(hvar1bin12wrong , "wrong classification", 'p')
	leg.Draw()

	#cvar1.SaveAs("test.png")
	cvar1.SaveAs("./plot/prune_test.pdf")

	dividers = cal_divider(ID3_tree, cvar1)
	for line in dividers:
		line.Draw()
	cvar1.SaveAs("./plot/prune_divide_test.pdf")

	print(hvar1bin12wrong.GetN())
	total_accuracy = 1 - 1.0*hvar1bin12wrong.GetN()/nentries
	print(total_accuracy)
	print('mistake 0 as 1: '),
	print(mistake0as1)
	print('mistake 1 as 0: '),
	print(mistake1as0)

	AUC = 1 - 0.5*(mistake0as1 + mistake0as1 + mistake1as0)/len(events)
	print('AUC: '),
	print(AUC)

	raw_input('Press Enter to continue...')

	return

main()

 


