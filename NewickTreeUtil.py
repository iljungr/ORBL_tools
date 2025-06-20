#!/usr/bin/env python
# Copyright 2025 Irwin Jungreis
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
Utilities for working with trees in Newick format (.nh).
Main class is NewickNode.
Includes parsing and writing to .nh files, traversal-related utilities,
and utilities for calculating probabilities under an evolutionary model.
"""
from __future__ import division
from __future__ import print_function
import sys, os
from copy import deepcopy

class NewickParseError(Exception) : pass

class NewickNode(object) :
    """A tree node with branch length to parent."""
    def __init__(self, name = '', distanceToParent = None) :
        self.children = []
        self.name = name
        self.distanceToParent = distanceToParent # None for root
    def is_leaf(self) :
        return self.children == []
    def aux_print_string(self) :
        """
        Override this in derived class, or set NewickNode.aux_print_string at runtime to
            to make print_tree print additional information for each node.
        """
        return None
    def copy(self, nodeClass = None, *initArgs) :
        """
        Return a deep copy of self. Args nodeClass and initArgs are obsolete but
        kept for compatibility.
        """
        return deepcopy(self)
    def __repr__(self) :
        """Although __repr__ is supposed to be unambiguous, it would be too
           dangerous to print the whole tree, as it could be very long."""
        return 'NewickNode(%s, %s, %d children)' % (self.name,
            '%g' % self.distanceToParent if self.distanceToParent is not None else None,
            len(self.children))
    def as_nh(self) :
        """Return a string containing the Newick respresentation of self."""
        return dump_nh_str(self, beautify = False).rstrip() # Remove final \n
    def __getitem__(self, key) :
        """Allow accessing children with node[index] or even node[:]"""
        return self.children.__getitem__(key)
    def __setitem__(self, key, value) :
        """Allow setting children with node[index] = """
        self.children.__setitem__(key, value)
    def equiv(self, other) :
        """
        Return True if other is equivalent to self, i.e.: has the same type, name, and
            distanceToParent, and has equiv children in same order.
        Don't require same value of the _descendants cache.
        WARNING: do not define __eq__ in this way because we use nodes as dictionary keys,
            which requires a __hash__ function consistent with __eq__, so __eq__ things
            would be treated as the same dictionary key, but we wouldn't want equiv things
            to be treated as the same dictionary key unless they are equal (for
            example, in a dictionary whose keys are all the subnodes of a tree, we might
            want different values for two leaves having no name and equal distanceToParent
            even though they are equiv).
        """
        if type(other) != type(self) : # Require same class
            return False
        if (self.name != other.name or self.distanceToParent != other.distanceToParent or
                len(self.children) != len(other.children)) :
            return False
        return all(c1.equiv(c2) for c1, c2 in zip(self.children, other.children))
    def iter_nodes(self, leafOnly = False) :
        """Iterate through descendants in pre-order including only leaves if requested."""
        if self.is_leaf() or not leafOnly :
            yield self
        for child in self :
            for descendant in child.iter_nodes(leafOnly) :
                yield descendant
    def iter_nodes_post(self, leafOnly = False) :
        """Iterate through descendants in post-order including only leaves if requested."""
        for child in self :
            for descendant in child.iter_nodes_post(leafOnly) :
                yield descendant
        if self.is_leaf() or not leafOnly :
            yield self
    def iter_parent_child_pairs(self, parentOfSelf = None) :
        """Iterate through descendants in pre-order, returning parent and child."""
        if parentOfSelf is not None :
            yield parentOfSelf, self
        for child in self :
            for pair in child.iter_parent_child_pairs(self) :
                yield pair
    def iter_leaves(self) :
        """Iterate through leaf descendants in depth-first order."""
        for descendant in self.iter_nodes(leafOnly = True) :
            yield descendant
    def leaf_names(self) :
        return leaf_names(self)
    def subtree_branch_length(self, names = None) :
        """Return the total branch length of the minimal subtree containing all
              nodes whose name is in names, or None if no nodes are in names.
           If names is None, use all leaves (i.e., total branch length of tree.)"""
        return subtree_branch_length(self, names)
    def relative_branch_length(self, names) :
        """Return relative branch length of minimal subtree containing names."""
        return relative_branch_length(self, names)
    def height(self, ignoreBranchLengths = False) :
        return calc_height(self, ignoreBranchLengths = ignoreBranchLengths)
    def prune(self, namesToKeep, copy = True) :
        return prune_tree(self, namesToKeep, copy)
    def dump_nh(self, fileName, beautify = True) :
        dump_nh(self, fileName, beautify )
    def set_descendants(self) :
        """
        Set node._descendants of self and each of its descendent nodes to the set of
            all leaf names under the node.
        Calling this first will make repeated calls to subtree_branch_length faster.
        WARNING: it is responsibility of caller to call this again when tree changes.
        """
        for node in self.iter_nodes_post() :
            if node.is_leaf() :
                node._descendants = {node.name}
            else :
                node._descendants = set.union(*[child._descendants for child in node])
    def print_tree(self, outFile = sys.stdout) :
        print_tree(self, outFile)
    @staticmethod
    def beautify_nh(inFileName, outFileName = None) :
        beautify_nh(inFileName, outFileName)

def print_tree(node, outFile = sys.stdout, indentLevel = 0) :
    """
    Print the tree in pre-order with indentation showing tree structure.
    By default, only distanceToParent and name are shown, but additional information
        can be shown by overriding aux_print_string in a derived class or setting
        NewickNode.aux_print_string at run time.
    Also, indicate which nodes are leaves (it can be deduced from tree structure, but
        it is easier to pick out visually if marked explicitly).
    """
    # Print the tree in pre-order
    distString = 'None' if node.distanceToParent is None else '%f' % node.distanceToParent
    nodeAuxString = node.aux_print_string()
    nodeDescription = '%s "%s"' % (distString, node.name)
    if nodeAuxString is not None :
        nodeDescription += ' "%s"' % nodeAuxString
    if node.is_leaf() :
        nodeDescription += ' Leaf'
    print('    ' * indentLevel + nodeDescription, file=outFile)
    for child in node.children :
        print_tree(child, outFile, indentLevel + 1)

def leaf_names(node) :
    # Return a list with the names of all leaves, in depth-first left-to-right order.
    return [leaf.name for leaf in node.iter_leaves()]

def calc_height(node, ignoreBranchLengths = False) :
    """
    Return height of the tree, i.e., length of longest path from a leaf up to the root.
    If ignoreBranchLengths, ignore distanceToParent and treat all edges as length 1.
    """
    if node.is_leaf() :
        return 0
    def one_level_height(child) :
        if ignoreBranchLengths :
            return 1
        elif child.distanceToParent is not None :
            return child.distanceToParent
        else :
            raise NewickParseError('distanceToParent is None. '
                                   'Use ignoreBranchLengths = True.')
    return max(calc_height(child, ignoreBranchLengths) + one_level_height(child)
               for child in node)

def calc_all_depths(root, ignoreBranchLengths = False) :
    """Return {node : depth, for all nodes in the tree}."""
    depthDict = {root : 0 if ignoreBranchLengths else 0.0}
    for node in root.iter_nodes() :
        nodeDepth = depthDict[node]
        for child in node.children :
            depthDict[child] = nodeDepth + (1 if ignoreBranchLengths else
                                            child.distanceToParent)
    return depthDict

def parse_nh_str(treeStr) :
    """
    Parse a string in Newick tree format and return the root.
    Although not part of the Newick specification, we treat lines at the beginning
        that start with '#' as comments, and ignore them.
    We do not handle Newick comments enclosed in square brackets.

    Newick format example:
    ((((((dmel:0.061361,(dsim:0.054894,dsec:0.031243):0.031837):0.063495,\
    (dyak:0.111338,dere:0.100461):0.039892):0.357431,dana:0.581114):0.243592,\
    (dpse:0.033045,dper:0.036095):0.495254):0.224541,dwil:0.801425):0.249420,\
    ((dvir:0.301255,dmoj:0.453117):0.141069,dgri:0.434875):0.249455);
    Names before colon for internal nodes are optional.
    """
    treeStr = treeStr.rstrip()
    while treeStr.startswith('#') : # Remove comment lines
        treeStr = treeStr[treeStr.find('\n') + 1 :]
    if treeStr[-1] != ';' :
        treeStr += ';' # Technically, it should have been there.
    treeStr = treeStr
    curToken = ''
    nodeStack = [NewickNode()]
    for ch in treeStr :
        if ch == '(' :
            nodeStack.append(NewickNode())
            nodeStack[-2].children.append(nodeStack[-1])
        elif ch in '),;' :
            nameAndDist = curToken.split(':')
            nodeStack[-1].name = nameAndDist[0].strip()
            if len(nameAndDist) > 1 :
                nodeStack[-1].distanceToParent = float(nameAndDist[1].strip())
            curToken = ''
            if ch == ';' :
                return nodeStack[0]
            if ch == ',' :
                nodeStack[-1] = NewickNode()
                nodeStack[-2].children.append(nodeStack[-1])
            else :
                del nodeStack[-1]
        else :
            curToken += ch

def parse_nh(infileOrName) :
    """Like parse_nh_str, but read tree from a file instead of a string."""
    try :
        with open(os.path.realpath(os.path.expanduser(infileOrName)), 'rt') as infile :
            treeStr = infile.read()
    except (TypeError, AttributeError) :  # It was a file, not a file name :
        treeStr = infileOrName.read()
    treeStr = treeStr.rstrip()
    return parse_nh_str(treeStr)

def dump_nh_str(root, beautify = True) :
    """Return a string containing the tree in Newick format."""
    def dump_nh_str_low(node, isRoot = False) :
        partialResult = ''
        if node.is_leaf() :
            partialResult += node.name
        else :
            partialResult += '('
            for childInd, child in enumerate(node.children) :
                partialResult += dump_nh_str_low(child)
                if childInd < len(node.children) - 1 :
                    partialResult += ','
            partialResult += ')'
            if node.name != '' :
                partialResult += node.name
        if not isRoot and node.distanceToParent is not None :
            partialResult += ':%r' % node.distanceToParent
        return partialResult
    result = dump_nh_str_low(root, isRoot = True) + ';\n'
    if beautify :
        result = beautify_nh_str(result)
    return result

def dump_nh(root, fileOrName, beautify = True) :
    """Write the tree starting at root into a file, in Newick format (.nh)."""
    resultStr = dump_nh_str(root, beautify)
    try :
        with open(os.path.realpath(os.path.expanduser(fileOrName)), 'wt') as outFile :
            outFile.write(resultStr)
    except (TypeError, AttributeError) : # It was a file, not a file name
        fileOrName.write(resultStr)

def beautify_nh_str(nhStr) :
    """Input is a string holding a tree in Newick format. Return a string with the same
       tree in Newick format with nice indentation. Preserve lines at the beginning that
       start with '#', if any."""
    lines = nhStr.split('\n')
    result = ''
    indent = 0
    for line in lines :
        if line.startswith('#') : # Pass comment lines through unchanged
            assert indent == 0, 'We only handle comments at the beginning.'
            result += line + '\n'
            continue
        for char in line :
            if char == '(' :
                indent += 1
                result += '(\n'
                result += '    ' * indent
            elif char == ')' :
                indent -= 1
                result += '\n'
                result += '    ' * indent
                result += ')'
            elif char == ',' :
                result += ',\n'
                result += '    ' * indent
            elif char not in ' \t' : # Ignore existing white space
                result += char
    result += '\n'
    return result

def beautify_nh(inFileOrName, outFileOrName = None) :
    """Read a tree in Newick format and write with nice indentation.
       If outFileName is None, overwrite inFileName (which must be a name, not a file)."""
    outFileOrName = inFileOrName if outFileOrName is None else outFileOrName
    try :
        with open(os.path.realpath(os.path.expanduser(inFileOrName)),  'rt') as inFile :
            inStr = inFile.read()
    except (TypeError, AttributeError) :  # It was a file, not a file name :
        assert outFileOrName != inFileOrName
        inStr = inFileOrName.read()
    result = beautify_nh_str(inStr)
    try :
        with open(os.path.realpath(os.path.expanduser(outFileOrName)), 'wt') as outFile :
            outFile.write(result)
    except (TypeError, AttributeError) :  # It was a file, not a file name :
        outFileOrName.write(result)

def prune_tree(tree, namesToKeep, copy = True) :
    """
    Return a tree containing the nodes with specified names and minimal set of ancestors.
       namesToKeep is any iterable.
       Combine consecutive branches when possible.
    If tree is a root (i.e., distanceToParent is None) then make new root also have
       distanceToParent None, otherwise make it the total distance from the new root to
       the parent of the original root.
    Return the new root, or None if there are no nodes left.
    If copy is false, change the original tree.
    """
    if not isinstance(namesToKeep, set) :
        # Convert iterator, list, etc to set for fast repeatable access
        namesToKeep = set(namesToKeep)
    if copy :
        tree = tree.copy()
    if tree.is_leaf() :
        return tree if tree.name in namesToKeep else None
    for childIndex, child in enumerate(tree.children) :
        tree.children[childIndex] = prune_tree(child, namesToKeep, copy = False)
    tree.children = [child for child in tree.children if child is not None]
    if tree.name not in namesToKeep and tree.children == [] :
        return None
    if len(tree.children) == 1 :
        if tree.distanceToParent is not None :
            tree.children[0].distanceToParent += tree.distanceToParent
        else :
            tree.children[0].distanceToParent = None
        tree = tree.children[0]
    return tree

def find_node_by_name(root, targetName) :
    """Return a node in the tree with targetName, or None."""
    for node in root.iter_nodes() :
        if node.name == targetName :
            return node
    return None

def path_to_node(root, nodeToFind) :
    """Return the indices in the children arrays starting at root that go to
       nodeToFind or None if nodeToFind is not in the tree."""
    if root is nodeToFind :
        return []
    if root.is_leaf() :
        return None
    for ind, child in enumerate(root) :
        pathFromChild = path_to_node(child, nodeToFind)
        if pathFromChild is not None :
            return [ind] + pathFromChild
    return None

def node_from_path(root, path) :
    """
    Given indices in the children arrays starting at root, return the node at which the
        path ends. Inverse of path_to_node.
    """
    node = root
    for ind in path :
        node = node.children[ind]
    return node

def common_ancestor(node, names = None, condition = None, leafOnly = True) :
    """
    Return lowest node under node that's an ancestor of all nodes with name in names
        or satisfying condition:node->bool, or None if there are no such nodes.
    Exactly one of names and condition must be None.
    If leafOnly, only check leaves.
    """
    assert (names is None) ^ (condition is None)
    if condition is None :
        condition = lambda n : n.name in names
    if node.is_leaf() :
        return node if condition(node) else None
    elif not leafOnly and condition(node) :
        return node
    subAncestors = [common_ancestor(child, condition = condition, leafOnly = leafOnly)
                    for child in node]
    subAncestors = [anc for anc in subAncestors if anc is not None]
    if len(subAncestors) == 0 :
        return None
    elif len(subAncestors) > 1 :
        return node
    else : # exactly 1
        return subAncestors[0]

def name_internal_nodes(root) :
    """Name every internal node of the tree by its pre-order ordinal number."""
    count = 0
    for node in root.iter_nodes() :
        if not node.is_leaf() :
            node.name = str(count)
            count += 1

def relative_branch_length(tree, names) :
    return subtree_branch_length(tree, names) / subtree_branch_length(tree)

def subtree_branch_length(tree, names = None) :
    """Return the total branch length of the minimal subtree containing all
          nodes whose name is in names.
       If names is None, use all leaves (i.e., total branch length of tree.)
       If calling repeatedly, call set_descendants first to make this faster."""
    leafNames = getattr(tree, '_descendants', None)
    if leafNames is None :
        leafNames = set(leaf_names(tree))
    if names is None :
        names = leafNames
    else :
        if any(name not in leafNames for name in names) :
            print('subtree_branch_length: Warning: names not in tree:',
                  [name for name in names if name not in leafNames], file=sys.stderr)
            # Remove extraneous names, especially ''
            names = [name for name in names if name in leafNames]
        names = set(names)
    bl = _subtree_branch_length_low(tree, names)[0]
    return 0.0 if bl is None else bl

def _subtree_branch_length_low(tree, names) :
    """Return a pair (branchLength, branchLengthWithRoot).
       branchLength is for minimal subtree including all names
       branchLengthWithRoot is for minimal subtree including all names and root"""
    if hasattr(tree, '_descendants') and len(names & tree._descendants) == 0 :
        return None, 0.0
    childBLpairs = [_subtree_branch_length_low(child, names) for child in tree.children]
    numNotNone = sum(childBL is not None for childBL, childBLwR in childBLpairs)
    if numNotNone == 0 :
        if tree.name in names :
            return 0.0, 0.0
        else :
            return None, 0.0
    branchLengthWithRoot = sum(childBLwR + child.distanceToParent
                               for child, (childBL, childBLwR) in zip(tree.children,
                                                                      childBLpairs)
                               if childBL is not None)
    if tree.name in names or numNotNone > 1 : # root is in the minimal tree
        branchLength = branchLengthWithRoot
    else :  # root is not in the minimal tree (sum below has only one element)
        branchLength = sum(childBL for childBL, childBLwR in childBLpairs
                           if childBL is not None)
    return branchLength, branchLengthWithRoot

def root_unrooted_tree(origRoot, childOfNewRoot = None) :
    """
    Given the root of a tree whose root has 3 children but is otherwise binary, and
        a node of that tree (childOfNewRoot), change the tree to a binary tree whose root
        is a new node on the branch going up from childOfNewRoot, and return the new root.
    Make the new root equidistant from its two children.
    If childOfNewRoot is None, use the result of find_best_pivot.
    """
    assert len(origRoot.children) == 3, 'Root does not have 3 children.'
    assert all(len(node.children) in [0, 2]
               for node in origRoot.iter_nodes()
               if node != origRoot), 'Rest of tree is not binary.'
    assert childOfNewRoot != origRoot

    if childOfNewRoot is None :
        childOfNewRoot = find_best_pivot(origRoot)

    parentDict = {childNode : parentNode
                  for parentNode, childNode in origRoot.iter_parent_child_pairs()}
    pivotNode = parentDict[childOfNewRoot]

    newRoot = NewickNode()
    newRoot.children.append(childOfNewRoot)
    newRoot.children.append(pivotNode)
    pivotNode.children.remove(childOfNewRoot)

    # Put newRoot halfway between its children
    childOfNewRoot.distanceToParent /= 2
    oldDist = pivotNode.distanceToParent # Save it for later
    pivotNode.distanceToParent = childOfNewRoot.distanceToParent

    while pivotNode != origRoot :
        nextPivot = parentDict[pivotNode]
        pivotNode.children.append(nextPivot)
        nextPivot.children.remove(pivotNode)
        nextPivot.distanceToParent, oldDist = oldDist, nextPivot.distanceToParent
        pivotNode = nextPivot

    return newRoot

def find_best_pivot(root) :
    """
    Given an unrooted tree, find the "best" rooting in the sense of minimizing the
    variance of the depths of the leaves. Return the node that will become the child
    of the new root.
    """
    pivotVariances = [] # [(varianceOfLeafDepths if node were childOfNewRoot, node)...]
    for node in root.iter_nodes() :
        if node == root :
            continue
        copiedTree = root.copy()
        copiedNode = node_from_path(copiedTree, path_to_node(root, node))
        newRoot = root_unrooted_tree(copiedTree, copiedNode)
        depthDict = calc_all_depths(newRoot)
        leafDepths = [d for n, d in depthDict.items() if n.is_leaf()]
        meanLeafDepth = sum(leafDepths) / len(leafDepths)
        varLeafDepth = sum((d - meanLeafDepth)**2 for d in leafDepths) / len(leafDepths)
        pivotVariances.append((varLeafDepth, node))
    return min(pivotVariances)[1]

def sort_tree(tree, namesInOrder, exact = True) :
    """
    namesInOrder is a list containing the leaf names of the tree or a subset thereof.
    Without changing the tree topology, reorder children as needed so that the depth-first
        order of the leaves whose names are in namesInOrder is as close as possible to
        their order in namesInOrder, and so that those leaves come before other leaves
        to the extent possible.
    Raise NameError if the set of namesInOrder includes names that are not in leaf names.
    If exact, raise ValueError if the order of names is incompatible with the tree
        topology (but tree will still be sorted to be as close as possible).
    NOTE: This calls set_descendants, which sets node._descendants for every node.
    """
    treeNames = tree.leaf_names()
    if (len(namesInOrder) != len(set(namesInOrder)) or
        len(treeNames)    != len(set(treeNames))    or
        len(set(namesInOrder) - set(treeNames)) != 0) :
        raise NameError
    nameIndexTable = {name : ii for ii, name in enumerate(namesInOrder)}
    tree.set_descendants()

    def ave_rank(aNode) :
        """
        Sort key for child nodes of some parent.
        Return the average rank in namesInOrder of the names of those descendants of aNode
            whose names are in namesInOrder; return infinity if there are none, so
            that nodes with names in namesInOrder will come before others when possible.
        """
        ranks = [nameIndexTable[name] for name in aNode._descendants
                 if name in nameIndexTable]
        return sum(ranks) / len(ranks) if ranks else float('inf')

    for node in tree.iter_nodes() :
        if not node.is_leaf() :
            node.children.sort(key = ave_rank)
    if exact :
        treeOrder = [name for name in tree.leaf_names()
                     if name in nameIndexTable]
        for name, request in zip(treeOrder, namesInOrder) :
            if name != request :
                msg = 'Order incompatible with tree topology. %s != %s' % (name, request)
                raise ValueError(msg)
