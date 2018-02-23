from abc import ABCMeta, abstractmethod
import utils
from collections import defaultdict, Counter
import numpy as np
import operator
from letterconfig import *


def init_back_table(dct):
    """Get back table for the current genetic code"""
    back_table = defaultdict(list)
    for aa, codon in zip(dct.forward_table.values(), dct.forward_table.keys()):
        back_table[aa].append(codon)
    return back_table


class AbsAncest:
    __metaclass__ = ABCMeta

    def __init__(self, tree, nodestates):
        self.tree = tree
        self.nodestates = nodestates
        self._set_tree()

    def _set_tree(self):
        nnode = 0
        for node in self.tree.traverse("postorder"):
            node.add_features(ind=nnode)
            nnode += 1
            # add nodestates infos to leaves
            if node.is_leaf():
                # only genome
                tmp_dct = self.nodestates.get(node.name, {})
                node.add_features(state=tmp_dct)
        # add total node count as a feature at the root
        self.tree.add_features(node_count=nnode)

    @abstractmethod
    def label_internal(self, state_mat, **kwargs):
        pass

    def mat_to_dict(self, charlist, ignore_this={}):
        """Return a dict from the mat"""
        mat_val = utils.makehash(1, list)
        for node in self.tree.traverse():
            for i, val in enumerate(charlist):
                for st in self.state_mat[i, node.ind]:
                    ignore_corr = (ignore_this.get(val, "") == st)
                    if not ignore_corr:
                        mat_val[node.name][st].append(val)
        return mat_val

    @classmethod
    def flip_rea_forward(clc, nodestates):
        """Flip rea data dict"""
        new_dt = utils.makehash(1, set)
        state_map = defaultdict(set)
        for (genome, aarea) in nodestates.items():
            nl = [(c, aa) for aa, codons in aarea.items() for c in codons]
            for (c, aa) in nl:
                new_dt[genome][c].add(aa)
                state_map[c].add(aa)

        return new_dt, state_map

    @classmethod
    def make_codonrea_matrice(clc, tree, nodestates, alphmap={}):
        new_dict, codon_map = clc.flip_rea_forward(nodestates)
        codon_list = sorted(codon_map.keys())
        # each entry will be a set
        # hmm not the best way, but hey I'm out of time
        tree_err = "Tree should be root and have total node count as feature"
        assert tree.is_root() and 'node_count' in tree.features, tree_err
        state_mat = np.zeros(
            (len(codon_list), tree.node_count), dtype='object')
        state_mat.fill(set([]))
        # looking only at leaf node
        for leaf in tree:
            leaf_state = new_dict[leaf.name]
            for (icod, state) in [(icod, set(leaf_state.get(cod, alphmap[cod]))) for icod, cod in enumerate(codon_list)]:
                state_mat[icod, leaf.ind] = state
        return state_mat, codon_map, codon_list


class LKLBasedAns(AbsAncest):

    class MLtype(object):
        marginal = 'margin'
        join = 'join'
        accepted_attr = [marginal, join]

    def __init__(self, tree, nodestates, mltype=MLtype.marginal):
        super(LKLBasedAns, self).__init__(tree, nodestates)
        if mltype not in LKLBasedAns.MLtype.accepted_attr:
            mltype = MLtype.marginal
        self.mltype = mltype

    def label_internal(self, state_mat, **kwargs):
        self.state_mat = state_mat
        if self.mltype == MLtype.marginal:
            return self._margin()
        else:
            return self._join()
        raise NotImplemented(
            "Likelihood based methods are not implemented yet")

    def _margin(self):
        print('margin')
        pass

    def _join(self):
        print('join')
        pass


class DolloParsimony(AbsAncest):

    def __init__(self, tree, nodestates, binary_mode=False, enforce_order=False, sort_by_size=False):
        super(DolloParsimony, self).__init__(tree, nodestates)
        self.binary_mode = binary_mode
        self.enforce_order = enforce_order
        self.sort_by_size = sort_by_size
        # only use enforce order if

    def label_internal(self, state_mat, header_list=None, header_map=None, null_map=None, **kwargs):
        if not self.binary_mode and None in (header_list, header_map):
            raise ValueError(
                "header_list and header_map are needed in non binary_mode")

        self.state_mat = state_mat
        self.header_list = header_list
        self.header_map = header_map
        self.null_map = null_map
        self._build_internal_label()
        return self.state_mat

    def __const_term_state(self, term_list, smat, char_states, nullstate):
        is_valid_for_c = {}
        if self.enforce_order:
            # here we suppose that
            char_states = reversed(sorted(list(char_states)))
        for c in char_states:
            # only consider positive state
            if c and c != nullstate:
                tmp_c = [sum([1 for x in smat[nlist] if c in x])
                         for nlist in term_list]
                if len([x for x in tmp_c if x > 0]) >= 2:
                    is_valid_for_c[c] = tmp_c

        # sort, whether by total number
        if self.enforce_order:
            return [(x, is_valid_for_c[x]) for x in char_states if x in is_valid_for_c.keys()]
        if self.sort_by_size:
            possible_order = sorted([(x, np.count_nonzero(y)) for x, y in is_valid_for_c.items(
            )], key=operator.itemgetter(1), reverse=True)
        else:
            # here we sort by subtree weight, ignoring excluded subgroup
            possible_order = sorted([(x, y[-1]) for x, y in is_valid_for_c.items()],
                                    key=operator.itemgetter(1), reverse=True)
        return possible_order

    def _build_internal_label(self):
        """Build internal state for dollo model
        see Farris, 1997 : http://sysbio.oxfordjournals.org/content/26/1/77.abstract
        """

        leaves_list_ind = set([l.ind for l in self.tree])
        get_poss_states = (lambda x: [0, 1]) if self.binary_mode else (
            lambda x: x[1][x[0]])

        def get_excluded_node_set(node):
            """ this function serve to get T-A according to Farris
            algorithm"""
            return list(leaves_list_ind - set([l.ind for l in node]))

        nullstate = 0
        for char_i, cur_char in enumerate(self.header_list):
            char_states = get_poss_states((cur_char, self.header_map))
            if not self.binary_mode:
                nullstate = self.null_map[cur_char]
            for node in self.tree.traverse("postorder"):
                # dollo look only at internal nodes
                if not node.is_leaf():
                    group_list = [[n.ind for n in child]
                                  for child in node.get_children()]
                    group_list.append(get_excluded_node_set(node))
                    pos_state = self.__const_term_state(group_list, self.state_mat[
                                                        char_i, :], char_states, nullstate)
                    # choose first state directly
                    # modify this to consider all potential solutions
                    chosen_state = {nullstate}
                    if pos_state:
                        chosen_state = {pos_state[0][0]}
                    # print chosen_state, node.name, node.ind, char_i, cur_char
                    self.state_mat[char_i, node.ind] = chosen_state


class FitchParsimony(AbsAncest):

    def __init__(self, tree, nodestates):
        super(FitchParsimony, self).__init__(tree, nodestates)

    def label_internal(self, state_mat, header_list=None, header_map=None, null_map=None, **kwargs):
        self.state_mat, self.header_list = state_mat, header_list
        # self.state_mat, self.codon_list = self.make_codonrea_matrice(self.tree, self.nodestates, self.alphmap)
        self._bottomup()
        self._updown()
        return self.state_mat

    def _bottomup(self):
        for i, st in enumerate(self.header_list):
            for node in self.tree.traverse("postorder"):
                if node.is_leaf():
                    # case previously resolved
                    pass
                else:
                    intersect = set.intersection(
                        *(self.state_mat[i, child.ind] for child in node.get_children()))
                    union = set.union(*(self.state_mat[i, child.ind]
                                        for child in node.get_children()))
                    if intersect:
                        self.state_mat[i, node.ind] = intersect
                    else:
                        self.state_mat[i, node.ind] = union

    def _get_node_optimal_state(self, node, ichar):
        accepted_state = self.state_mat[ichar, node.ind]
        char_score = {}
        for c in accepted_state:
            in_subtree = sum(
                [(1 if c in self.state_mat[ichar, n.ind] else 0) for n in node])
            out_subtree = sum([(1 if c in self.state_mat[ichar, n.ind] else 0)
                               for n in node.get_tree_root() if n not in node])
            char_score[c] = in_subtree - out_subtree
        return {Counter(char_score).most_common()[0][0]}

    def _updown(self):
        for ichar, st in enumerate(self.header_list):
            # choose state for root
            if len(self.state_mat[ichar, self.tree.ind]) > 1:
                rstate = self._get_node_optimal_state(self.tree, ichar)
                self.state_mat[ichar, self.tree.ind] = rstate

            for node in self.tree.iter_descendants("preorder"):
                parent_state = self.state_mat[ichar, node.up.ind]
                if parent_state.intersection(self.state_mat[ichar, node.ind]):
                    self.state_mat[ichar, node.ind] = parent_state
                else:
                    self.state_mat[ichar, node.ind] = self._get_node_optimal_state(
                        node, ichar)


class SingleNaiveRec(object):
    """A NaiveFitch algorithm for finding the most parcimonious solution"""

    def __init__(self, tree, reassigned, ori_aa, dest_aa, dct, codon_rea=(None, None), mode="fitch"):
        self.id = {}
        self.tree = tree
        self.corr = {'0': ori_aa, '1': dest_aa}
        self.codon_rea_global, self.codon_rea_filtered = codon_rea
        colors = ['#a6cee3', '#1f78b4', '#b2df8a',
                  '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f']
        self.dct = dct
        self.back_table = init_back_table(dct)
        codon_list = self.back_table[aa_letters_3to1[ori_aa]]
        self.colors = dict(zip(codon_list, colors))
        for leaf in tree:
            if leaf.name in reassigned:
                leaf.add_features(reassigned={1})
                leaf.add_features(rea=self.corr['1'])
                leaf.add_features(state=dest_aa)
            else:
                leaf.add_features(reassigned={0})
                leaf.add_features(rea=self.corr['0'])
                leaf.add_features(state=ori_aa)

        self.ori_aa = ori_aa
        self.dest_aa = dest_aa
        self.ori_aa1 = aa_letters_3to1[ori_aa]
        self.dest_aa1 = aa_letters_3to1[dest_aa]
        self.newick = tree.write(features=['name', 'dist', 'support', 'state'])
        self.mode = mode
        self._bottomup(self.tree)
        self._topdown(self.tree)

    def update_codon_data(self, codon_rea):
        """Update codon reassignment data (global and filtered)
        """
        self.codon_rea_global, self.codon_rea_filtered = codon_rea

    def write_tree(self, outfile):
        """Export newick to  a file"""
        with open(outfile, 'w') as OUT:
            OUT.write(self.newick)

    def is_valid(self, thresh=1):
        """Naive function to test whether or not if the current reassignment is valid"""
        tmptree = self.tree.copy()
        for l in tmptree.traverse():
            if not l.is_leaf():
                l.del_feature('reassigned')
                l.del_feature('rea')
            elif (l.state == self.dest_aa and ('lost' is l.features and l.lost)) or ('count' in l.features and l.count < thresh):
                l.add_features(reassigned={0})
                l.add_features(rea=self.corr['0'])
                l.add_features(state=self.ori_aa)
            elif ('lost' in l.features and not l.lost) and l.state == self.ori_aa and ('count' in l.features and l.count >= thresh):
                l.add_features(reassigned={1})
                l.add_features(rea=self.corr['1'])
                l.add_features(state=self.dest_aa)
        # At least one node with reassignment persist in the data
        for node in tmptree:
            if node.rea == self.dest_aa:
                return True
        return False

    @classmethod
    def _fitch(clc, tree, corr):
        """Fitch algorithm part 1 : bottom up"""
        for node in tree.traverse('postorder'):
            if not node.is_leaf():
                intersect = []
                union = set()
                for c in node.get_children():
                    union = union | c.reassigned
                    intersect.append(c.reassigned)
                intersect = set.intersection(*intersect)
                if(intersect):
                    node.add_features(reassigned=intersect)
                else:
                    node.add_features(reassigned=union)

                node.add_features(
                    rea="/".join([corr[str(r)] for r in node.reassigned]))

    @classmethod
    def _dollo(clc, tree):
        nodelist = tree.get_leaves()
        for node in tree.traverse('postorder'):
            if not node.is_leaf():
                internal_grp = [max([next(iter(c.reassigned))
                                     for c in n]) for n in node.get_children()]
                if not node.is_root():
                    internal_grp.append(max(
                        [next(iter(c.reassigned)) for c in nodelist if c.name not in node.get_leaf_names()]))
                if np.count_nonzero(internal_grp) >= 2:
                    node.add_features(reassigned={1})
                else:
                    node.add_features(reassigned={0})

    def _bottomup(self, tree):
        if self.mode == "fitch":
            self._fitch(tree, self.corr)
        elif self.mode == "dollo":
            self._dollo(tree)
        else:
            raise NotImplementedError(
                "The method %s you asked for is not implemented" % self.mode)

    def _topdown(self, tree):
        """Fitch algorithm part 2 : top down"""
        # Since it's not need, just pass
        pass

    def is_reassigned(cls, node, strict=True):
        """ return True if a node has undergoned reassignment """
        if "reassigned" in node.features and (node.reassigned == {1} or (not strict and 1 in node.reassigned)):
            return True
        return False

    def get_species_list(self, limit_to_suspected=False):
        """Get the current species list for this tree"""
        slist = set()
        if limit_to_suspected:
            for node in self.tree.traverse():
                if 'reassigned' in node.features and node.reassigned == {1}:
                    slist.update(node.get_leaf_names())
        else:
            slist = set(self.tree.get_tree_root().get_leaf_names())
        return slist

    def get_distance_to_rea_node(self, node):
        """Get the distance of node to the closest reassigned lca"""
        cur_node = node
        if isinstance(node, str):
            cur_node = self.tree & node
        dist = 0
        while cur_node is not None and not self.is_reassigned(cur_node, strict=False):
            dist += 1
            cur_node = cur_node.up
        return dist

    def has_codon_data(self):
        # base codon_data on filtered position only
        return self.codon_rea_filtered is not None
