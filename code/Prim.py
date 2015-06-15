#!/bin/python3 -tt
################################################################################
# Name:			Peter Gordon
# Email:		peter.gordon@csu.fullerton.edu
# Course:		MATH 370 (Mathematical Model Building)
#				Mon./Weds. 2:30--3:45 PM
# Instructor:	Dr. L. Smith
################################################################################
# (C) 2014 Peter Gordon <peter.gordon@csu.fullerton.edu>
# This source code is released under the terms of the Apache License 2.0;
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# These data structures and related functions/methods are based heavily on
# the csuf::Edge, csuf::Graph, and csuf::prim C++ code created by Dr. Kevin
# Wortman, also under the Apache 2.0 license, and available online at
#	http://mirror.thecodergeek.com/src/graph.h
#	http://mirror.thecodergeek.com/src/prim.h
################################################################################
# This module implements Prim's MST algorithm.
################################################################################

from heapq import heappop, heappush
from Graph import Edge, Graph, Vertex
from operator import itemgetter


def bridgeEnd(theGraph, inTree, edgeID):
	vert1 = theGraph.getEdge(edgeID).fromVert
	vert2 = theGraph.getEdge(edgeID).toVert
	if not inTree[vert1]:
		vert1, vert2 = vert2, vert1
	if not (inTree[vert1] and not inTree[vert2]):
		return False
	return (vert1, vert2)


def visitVertex(heap, inTree, theGraph, vertID):
	inTree[vertID] = True
	print("DEBUG visitVertex: @ vertex {}".format(theGraph._verts[vertID].name))
	for edgeID in theGraph.getIncidentEdges(vertID):
		if bridgeEnd(theGraph, inTree, edgeID):
			print("DEBUG visitVertex: found a bridge end!")
			if edgeID not in map(itemgetter(1), heap):
				print("DEBUG visitVertex: Adding edge to PQ: {e}".format(e=theGraph.edgeStr(edgeID)))
				heappush(heap, (theGraph.getEdge(edgeID).weight, edgeID))


def doPrimMST(theGraph):
	n = theGraph.numVerts()
	assert n > 0, "Empty graph."

	heap = []
	inTree = [False] * n
	MSTedges = []
	visitVertex(heap, inTree, theGraph, 0)
	while len(MSTedges) < (n-1):
		alreadySeen = False
		while not alreadySeen:
			(ignored, edgeID) = heappop(heap)
			print("DEBUG Prim: Popped edge from PQ: {e}".format(e=theGraph.edgeStr(edgeID)))
			alreadySeen = bridgeEnd(theGraph, inTree, edgeID)
			if alreadySeen:
				inVert, outVert = alreadySeen
		visitVertex(heap, inTree, theGraph, outVert)
		print("DEBUG Prim: Add edge to MST: {e}".format(e=theGraph.edgeStr(edgeID)))
		MSTedges.append(edgeID)

	return theGraph.getSubgraph(MSTedges)

