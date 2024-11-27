from dataclasses import dataclass
from typing import List, Any, Callable, Dict, Tuple
import random
import numpy as np

DiscreteVector = List[int]

@dataclass
class BKTreeNode:
    # The vector that is compared via the distance metric
    vector: DiscreteVector
    # The list of elements that are associated with the above vector
    elements: List[Any]
    # An id associated with the above vector
    id: int
    # The map from distance (e.g., D=1) to the child node associated with that distance
    children: Dict[int, "BKTreeNode"]
    
    @staticmethod
    def make_empty(vector=[]) -> "BKTreeNode":
        return BKTreeNode(vector=vector, elements=[], id = -1, children={})

    def is_empty(self) -> bool:
        return not self.vector

# This turned out to be not at all helpful -- just building the BKTree from a random start node is fine
def find_best_start(vectors: List[DiscreteVector],
                    distance: Callable[[DiscreteVector, DiscreteVector], int]) -> DiscreteVector:
    """
    Find the "most central" vector to be used for constructing a BK-tree's root node.
    """
    M = 10
    indices = [random.randint(0, len(vectors)-1) for _ in range(M)]
    best_stddev = 0.0
    best_vector = vectors[0]
    for row in vectors:
        dists = []
        for i in indices:
            dists.append(distance(vectors[i], row))
        stddev = np.std(dists)
        if stddev > best_stddev:
            best_stddev = stddev
            best_vector = row
    return best_vector

def bk_tree_insert(root_node: BKTreeNode,
                   elements: List[Any],
                   vector: DiscreteVector,
                   distance: Callable[[DiscreteVector, DiscreteVector], int]) -> BKTreeNode:
    """
    Insert the element into the BK tree rooted at root_node.

    :param root_node: The root of the BK tree.
    :param element: The element (ID, object containing data, whatever...) to be inserted.
    :param vector: The DiscreteVector associated with the element.
    :param distance: The distance callback for the DiscreteVector (e.g., hamming distance implementation).
    :returns: The node which contains element after insertion.
    """
    if root_node.is_empty():
        root_node.vector = vector
        root_node.elements.extend(elements)
        return root_node
    cur_node = root_node
    while cur_node is not None:
        k = distance(cur_node.vector, vector)
        if k == 0:
            cur_node.elements.extend(elements)
            return cur_node
        new_node = cur_node.children.get(k)
        if new_node is None:
            new_node = BKTreeNode(vector=vector, elements=list(elements), children={})
            cur_node.children[k] = new_node
            return new_node
        cur_node = new_node

def bk_tree_insert2(root_node: BKTreeNode,
                   elements: List[Any],
                   vector: DiscreteVector,
                   distance: Callable[[DiscreteVector, DiscreteVector], int]) -> BKTreeNode:
    """
    Insert the element into the BK tree rooted at root_node.

    :param root_node: The root of the BK tree.
    :param element: The element (ID, object containing data, whatever...) to be inserted.
    :param vector: The DiscreteVector associated with the element.
    :param distance: The distance callback for the DiscreteVector (e.g., hamming distance implementation).
    :returns: The node which contains element after insertion.
    """
    best_dist = 2**32
    if root_node.is_empty():
        root_node.vector = vector
        root_node.elements.extend(elements)
        return best_dist
    cur_node = root_node
    while cur_node is not None:
        k = distance(cur_node.vector, vector)
        best_dist = min(best_dist,k)
        if k == 0:
            cur_node.elements.extend(elements)
            return best_dist
        new_node = cur_node.children.get(k)
        if new_node is None:
            new_node = BKTreeNode(vector=vector, elements=list(elements), children={})
            cur_node.children[k] = new_node
            return best_dist
        cur_node = new_node


def bk_tree_lookup(root_node: BKTreeNode,
                   vector: DiscreteVector,
                   distance: Callable[[DiscreteVector, DiscreteVector], int],
                   skip: List[Any]) -> Tuple[List[BKTreeNode], int]:
    """
    Lookup the nearest neighbor(s) to the given vector and return their nodes and the best distance

    :param root_node: The root of the BK-tree.
    :param vector: The query DiscreteVector.
    :param skip: If a node only has the elements in this list then we skip this node
    :param collect_all: If True, find the distance of the nearest neighbor and then return all
        nodes that have that distance from the query vector. If False, just return the first one
        encountered.
    """
    if root_node.is_empty():
        return set()
    node_list = [root_node]
    results = []
    dist_best = 2**32
    while node_list:
        node = node_list.pop()
        dist = distance(node.vector, vector)

        # Check to see if the node contains only elements we want to skip (contains only close relatives)
        has_only_relatives = True
        for element_to_check in node.elements:
            if not element_to_check in skip:
                has_only_relatives = False
                break

        # We do not return nodes that have no elements or if the node only contains elements we want to skip
        ignore_this_node = (len(node.elements) == 0) or (has_only_relatives)
            
        # We return if we find an exact match (that isn't something we want to skip)
        if (dist == 0) and (not ignore_this_node):
            results = [node]
            dist_best = 0
            return results, dist_best
            #print(exact_matches)
        
        if not ignore_this_node:
            if dist < dist_best:
                # We've got a strictly better distance, clear the old set
                results = [node]
                dist_best = dist
            elif dist == dist_best:
                results.append(node)
        for next_dist, next_node in node.children.items():
            # Everything at or below next_node has dist(node.vector, ?) = next_dist.
            bound = abs(next_dist - dist)
            if (bound <= dist_best):
                node_list.append(next_node)
    return results, dist_best