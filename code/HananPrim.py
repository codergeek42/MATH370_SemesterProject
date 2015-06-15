#!/bin/python3 -tt
################################################################################
# Author:		Peter Gordon
# Email:		peter.gordon@csu.fullerton.edu
# Course:		MATH 370 (Mathematical Model Building)
#				Mon./Weds. 2:30--3:45 PM
# Instructor:	Dr. L. Smith
# Group with:	Casey Cao-Son and Brennan Truong
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
#################################################################################
# This implements the Hanan/Prim Hybridization approach to find a
# 1.5-approximation (with respect to total weight) of the MRST of a given set of
# points.
#################################################################################

from Graph import Edge, Graph, Vertex, makeGraphComplete, L1Dist
from itertools import product, tee
from Prim import visitVertex, bridgeEnd, doPrimMST
from heapq import heappop, heappush, heapify, merge
from operator import itemgetter
from math import pow

import matplotlib.pyplot as plot

class NaiveHananSolver:
	__slots__ = ("_graph", "_startingVerts", "_startingVertCoords", "_stationWeightFunc")

	def __init__(self, startingVerts):
		self._graph = Graph()
		Xs = Ys = set()

		self._startingVerts = []
		self._stationWeightFunc = lambda g,v: 0
		self._startingVertCoords = startingPts = [(vert.X, vert.Y) for vert in startingVerts]
		for vert in startingVerts:
			self._startingVerts.append(self._graph.addVertex(Vertex(vert.X, vert.Y, vert.name)))
			Xs = Xs.union([vert.X])
			Ys = Ys.union([vert.Y])
		hananPts = product(*(zip(*startingPts)))
		for (hananID, (vertX, vertY)) in enumerate(set(hananPts).difference(startingPts)):
			self._graph.addVertex(Vertex(vertX, vertY, "h"+str(hananID)))

		## NB: Make sure these all are in-scope, as we'll need their end-of-loop values below.
		x = rightX = y = upY = thisVertID = rightVertID = upVertID = rightUpVertID = None

		Xs = sorted(Xs)
		Ys = sorted(Ys)
		XsCurr, XsRight = tee(Xs)
		next(XsRight, None)
		for (x, rightX) in zip(XsCurr, XsRight):
			YsCurr, YsUp = tee(Ys)
			next(YsUp, None)
			for (y, upY) in zip(YsCurr, YsUp):
				thisVertID = self._graph.getVertexID(x, y)
				rightVertID = self._graph.getVertexID(rightX, y)
				upVertID = self._graph.getVertexID(x, upY)
				self._graph.addEdgeWithVertsAndWeight(thisVertID, rightVertID, abs(rightX-x))
				self._graph.addEdgeWithVertsAndWeight(thisVertID, upVertID, abs(upY-y))
				rightUpVertID = self._graph.getVertexID(rightX, upY)
				self._graph.addEdgeWithVertsAndWeight(rightVertID, rightUpVertID, abs(upY-y))
			self._graph.addEdgeWithVertsAndWeight(upVertID, rightUpVertID, abs(rightX-x))



	## stationWeightFunc should accept the graph object as the first parameter and the vertex ID
	## as the second.
	def addStationWeightFunc(self, stationWeightFunc):
		self._stationWeightFunc = stationWeightFunc


	def solve(self):
		def _areVerticesColinear(theGraph, vertIDs):
			assert len(vertIDs) == 3, "Expected 3 vertex IDs."

			(vert1, vert2, vert3) = (theGraph.getVertex(vertID) for vertID in vertIDs)
			if vert1.X == vert2.X == vert3.X:
				return "Y"
			if vert1.Y == vert2.Y == vert3.Y:
				return "X"
			return False


		def _removePassthroughPoints(theGraph, distFunc):
			## Each pass will remove all 2-order points that are not starting vertices,
			## where the points are essentially unneeded. For example, if there are
			## edges A-->B and B-->C, and B is not a starting vertex and has no other
			## neighbors, and A, B, and C are all colinear, replaces both of these
			## edges with a single edge A-->C and removes B from the graph.
			nothingRemoved = False
			while not nothingRemoved:
				nothingRemoved = True
				for (vertID, vertObj) in theGraph.getVerts():
					if theGraph.getDegreeOfVertex(vertID) == 2:
						nbrs = theGraph.getNeighborhood(vertID)
						if (_areVerticesColinear(theGraph, [vertID] + nbrs)
								and vertID not in self._startingVerts):
							(otherVert1ID, otherVert1Obj), (otherVert2ID, otherVert2Obj) = [
									(vertID, theGraph.getVertex(vertID)) for vertID in nbrs]
							theGraph.removeVertex(vertID)
							theGraph.addEdgeWithVertsAndWeight(otherVert1ID, otherVert2ID, 
									distFunc(otherVert1Obj, otherVert2Obj))
							nothingRemoved = False

		def _removeUselessPoints(theGraph, maxOrder, removeNearStartingVerts=False):
			## Remove all n-order or less Hanan points until grid unchanged.
			nothingRemoved = False
			for vertID in self._startingVerts:
				print ("DEBUG: solve: startingVerts is: ({ID}) {name}".format(ID=vertID, name=theGraph._verts[vertID].name))

			while not nothingRemoved:
				print("DEBUG: solve: Starting a {n}-order removal pass.".format(n=maxOrder))

				nothingRemoved = True
				for (vertID, vertObj) in theGraph.getVerts():
					print("DEBUG: solve: deg = {deg} @vertex {name}".format(name=vertObj.name,
						deg=theGraph.getDegreeOfVertex(vertID)))
					isAdjStartVert = len(set(theGraph.getNeighborhood(vertID)).intersection(self._startingVerts)) > 0
					if (theGraph.getDegreeOfVertex(vertID) <= maxOrder
							and vertID not in self._startingVerts
							and (isAdjStartVert if removeNearStartingVerts else True)):
						print("DEBUG: solve: @vertex ({ID}) {name}".format(ID=vertID, name=theGraph._verts[vertID].name))
						theGraph.removeVertex(vertID)
						nothingRemoved = False



		## This is the Naive approach: Simply take an MST of the Hanan Grid and remove
		## branches that don't terminate in a starting vertex.

		steinerTree = doPrimMST(self._graph)
		_removeUselessPoints(steinerTree, 1, False)
		_removePassthroughPoints(steinerTree, L1Dist)
		self._graph = steinerTree
		steinerTree.printGraph()

	def plotGraph(self):
		graphWeight = self.totalWeight()
		numStartingVerts = len(self._startingVerts)
		for (steinerPtID, steinerVert) in self._graph._verts.items():
			steinerPlot = plot.scatter(x=steinerVert.X, y=steinerVert.Y, c=(1,0,0), marker="o",
					linewidths=(4,))
			plot.annotate(s=steinerVert.name, xy=(steinerVert.X, steinerVert.Y), xytext=(-10,-10),
					textcoords="offset points", ha="right")

		for (startX, startY) in self._startingVertCoords:
			startingPlot = plot.scatter(x=startX, y=startY, c=(0,0,0), marker="s", linewidths=(8,))


		for (edgeID, edgeObj) in self._graph.getEdges():
			srcVert = self._graph.getVertex(edgeObj.fromVert)
			destVert = self._graph.getVertex(edgeObj.toVert)
			plot.plot((srcVert.X, destVert.X), (srcVert.Y, destVert.Y), "k-")


		plot.legend((startingPlot, steinerPlot), ("Starting Vertices", "Steiner Points"))
		plot.show()


	def totalWeight(self):
		return self.graphStats(printOutput=False)


	def graphStats(self, printOutput=True):
		if printOutput:
			print("Graph Vertices (total {n}):".format(n=self._graph.numVerts()))
		totalWeight = 0
		for (vertID, vertObj) in self._graph._verts.items():
			if printOutput:
				print("\t({ID}) {vertInfo!s}".format(ID=vertID, vertInfo=vertObj))
			if (vertObj.X, vertObj.Y) not in self._startingVertCoords:
				totalWeight += (self._stationWeightFunc)(self._graph, vertID)

		if printOutput:
			print("Graph Edges (total {m})".format(m=self._graph.numEdges()))
		for (edgeID, edgeObj) in self._graph._edges.items():
			if printOutput:
				print("\t({ID}) {src!s} -- {dest!s}  [{weight:.2f}]".format(
						ID=edgeID, src=self._graph._verts[edgeObj.fromVert],
						dest=self._graph._verts[edgeObj.toVert], weight=edgeObj.weight))
			totalWeight += edgeObj.weight

		if printOutput:
			print("Total weight: {weight:.2f}".format(weight=totalWeight))
		return totalWeight


if __name__ == "__main__":
	a = Vertex(0, 15, "A")
	b = Vertex(5, 20, "B")
	c = Vertex(16,  24, "C")
	d = Vertex(20, 20, "D")
	e = Vertex(33, 25, "E")
	f = Vertex(23, 11, "F")
	g = Vertex(35, 7, "G")
	h = Vertex(25, 0, "H")
	i = Vertex(10, 3, "I")
	grid = NaiveHananSolver([a,b,c,d,e,f,g,h,i])
	grid.addStationWeightFunc(lambda g,v: 1.2 * pow(g.getDegreeOfVertex(v), 1.5))
	grid.solve()
	grid.graphStats(printOutput=True)
	grid.plotGraph()

