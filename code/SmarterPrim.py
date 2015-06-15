#!/bin/python3 -tt
################################################################################
# Author:		Peter Gordon
# Email:		peter.gordon@csu.fullerton.edu
# Course:		MATH 370 (Mathematical Model Building)
#				Mon./Weds. 2:30--3:45 PM
# Instructor:	Dr. L. Smith
# Group:		Casey Cao-Son and Brennan Truong
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
# This implements the Modified Prim's Algorithm approach to find a
# 1.5-approximation (with respect to total weight) of the MRST of a given set of
# points.
#################################################################################

from Graph import Edge, Graph, Vertex, makeGraphComplete, L1Dist
from itertools import product, tee
from Prim import visitVertex, bridgeEnd, doPrimMST
from heapq import heappop, heappush, heapify, merge
from operator import itemgetter, attrgetter
from math import pow
import matplotlib.pyplot as plot


class SmarterPrimSolver:
	__slots__ = (
		"_graph", 				# The internal Graph structure
		"_startingVerts", 		# The list of starting vertices
		"_startingVertCoords",  # The coordinate pairs of the starting vertices
		"_rangeX", "_rangeY", 	# The range of values for the X and Y coordinates
		"_centerVert",			# A virtual central Vertex. Not actually a part
								# of the graph but  useful to heuristically
								# choose the better of two possible Steiner points.
		"_stationWeightFunc")	# The weight of each station

	def __init__(self, startingVerts):
		self._graph = Graph()
		Xs = Ys = set()

		self._stationWeightFunc = lambda g,v: 0
		self._startingVertCoords = []
		self._startingVerts = []
		startingPts = [(vert.X, vert.Y) for vert in startingVerts]

		for vert in startingVerts:
			self._startingVerts.append((self._graph.addVertex(Vertex(vert.X, vert.Y, vert.name)), vert.name))
			self._startingVertCoords.append((vert.X, vert.Y))
			Xs = Xs.union([vert.X])
			Ys = Ys.union([vert.Y])
		# hananPts = product(*(zip(*startingPts)))
		# for (hananID, (vertX, vertY)) in enumerate(set(hananPts).difference(startingPts)):
		#	self._graph.addVertex(Vertex(vertX, vertY, "h"+str(hananID)))
		#
		## NB: Make sure these all are in-scope, as we'll need their end-of-loop values below.
		# x = rightX = y = upY = thisVertID = rightVertID = upVertID = rightUpVertID = None

		self._rangeX = Xs = sorted(Xs)
		self._rangeY = Ys = sorted(Ys)
		self._centerVert = Vertex(sum(Xs)/float(len(Xs)), sum(Ys)/float(len(Ys)), "<CENTER>")


	## stationWeightFunc should accept the graph object as the first parameter and the vertex ID
	## as the second. For example
	def addStationWeightFunc(self, stationWeightFunc):
		self._stationWeightFunc = stationWeightFunc


	def solve(self):
		def _updateDistances(steinerTree, disconnectedVertDistDict, distFunc):
			for graphVertID in disconnectedVertDistDict:
				tempDists = []
				for (steinerVertID, steinerVertObj) in steinerTree.getVerts():
					dist = distFunc(steinerVertObj, self._graph.getVertex(graphVertID))
					tempDists.append((dist, steinerVertID))
				disconnectedVertDistDict[graphVertID] = min(tempDists)


		def _areVerticesColinear(steinerTree, vertIDs):
			assert len(vertIDs) == 3, "Expected 3 vertex IDs."

			(vert1, vert2, vert3) = [steinerTree.getVertex(vertID) for vertID in vertIDs]
			if vert1.X == vert2.X == vert3.X:
				return "Y"
			if vert1.Y == vert2.Y == vert3.Y:
				return "X"
			return False


		## Remove overlapping edges, to be run after the modified Prim MST algorithm.
		## Compensates for overcalculated edges in that loop by replacing them with shorter
		## indirect paths. For example, if the graph has colinear points
		## A--B--C, and edges A-->C and B-->C, replace the A-->C edge with A-->B.
		## (Otherwise the weight for B-->C is counted twice in the final total cost.)
		def _removeOverlaps(steinerTree, distFunc):
			nothingChanged = False
			while not nothingChanged:
				nothingChanged = True
				for (vertID, vertObj) in steinerTree.getVerts():
					nbrVerts = [(nbrVertID, steinerTree.getVertex(nbrVertID)) for nbrVertID in steinerTree.getNeighborhood(vertID)]
					nbrVertsByDist = [(nbrVertID, distFunc(nbrVertObj, vertObj)) for (nbrVertID, nbrVertObj) in nbrVerts]
					for (nbrVertID, ignored) in sorted(nbrVertsByDist, key=itemgetter(1)):
						for nextNbrVertID in steinerTree.getNeighborhood(nbrVertID):
							if nextNbrVertID == vertID:
								continue
							colinear = _areVerticesColinear(steinerTree, [vertID, nbrVertID, nextNbrVertID])
							verts = [(vertID, vertObj)] + [(tmpVertID, steinerTree.getVertex(tmpVertID))
									for tmpVertID in [nbrVertID, nextNbrVertID]]
							if not colinear:
								continue
							coord = attrgetter(colinear)
							verts.sort(key=lambda k: coord(k[1]))

							(vert1ID, vert1Obj) = verts[0]
							(vert2ID, vert2Obj) = verts[1]
							(vert3ID, vert3Obj) = verts[2]
							if ((vert1ID, vert1Obj.name) in self._startingVerts
									or (vert2ID, vert2Obj.name) in self._startingVerts
									or (vert3ID, vert3Obj.name) in self._startingVerts):
								continue
							wholeEdge = steinerTree.getEdgeBetween(vert1ID, vert3ID)
							if wholeEdge:
								steinerTree.removeEdge(wholeEdge.edgeID)
								steinerTree.addEdgeWithVertsAndWeight(vert1ID, vert2ID, distFunc(vert1Obj, vert2Obj))
								steinerTree.addEdgeWithVertsAndWeight(vert2ID, vert3ID, distFunc(vert2Obj, vert3Obj))
								nothingChanged = False


		## Adds the given vertex to the growing Steiner Tree by connecting it to the given
		## bridge vertex, inserting a Hanan point if needed.
		def _connectVertex(steinerTree, newVert, bridgeVertID, weightFunc):
			bridgeVert = steinerTree.getVertex(bridgeVertID)
			newVertID = steinerTree.addVertex(newVert)
			if (bridgeVert.X != newVert.X and bridgeVert.Y != newVert.Y):
				newSteinerVert = Vertex(bridgeVert.X, newVert.Y, "s"+str(steinerTree.numVerts()))
				newSteinerVert2 = Vertex(newVert.X, bridgeVert.Y, "s"+str(steinerTree.numVerts()))
				if weightFunc(newSteinerVert, self._centerVert) > weightFunc(newSteinerVert2, self._centerVert):
					newSteinerVert = newSteinerVert2
				newBridgeVertID = steinerTree.addVertex(newSteinerVert)
				steinerTree.addEdgeWithVertsAndWeight(newBridgeVertID, bridgeVertID, weightFunc(bridgeVert, newVert))
				bridgeVert = newSteinerVert
				bridgeVertID = newBridgeVertID
			steinerTree.addEdgeWithVertsAndWeight(bridgeVertID, newVertID, weightFunc(bridgeVert, newVert))



		## This is the Modified Prim's Algorithm, creating Hanan points on-demand.
		## NOTE: Vertex IDs in the Steiner Tree might differ from those of the same points in the
		## augmented original graph!
		n = self._graph.numVerts()
		assert n > 1, "Empty or single-element graph."


		allVerts = self._graph.getVerts()
		((ignored, firstVert), allVerts) = allVerts[0], allVerts[1:]
		steinerTree = Graph()
		steinerTree.addVertex(Vertex(firstVert.X, firstVert.Y, firstVert.name))

		## Continually tracks the minimum distance of each remaining point to the growing tree.
		## as a map of vertex ID (in graph) to a tuple of (distance, nearest steiner vertex ID)
		disconnectedDists = dict((vertID, (0,None)) for (vertID, vertObj) in allVerts)


		while len(disconnectedDists) > 0:
			print("DEBUG: solve: Finding next closest point...")
			_updateDistances(steinerTree, disconnectedDists, L1Dist)
			for (vertID, (dist, bridgeVertID)) in disconnectedDists.items():
				print("DEBUG: solve: Candidate is {f} -- {t} [{d}]".format(
						f=self._graph.getVertex(vertID), t=steinerTree.getVertex(bridgeVertID),
						d=dist))
			nextVertID = min(disconnectedDists, key=lambda idx: disconnectedDists[idx][0])
			(ignored, bridgeVertID) = disconnectedDists.pop(nextVertID)
			print("DEBUG: solve: Connecting {v1} to {v2}".format(
					v1=steinerTree.getVertex(bridgeVertID), v2=self._graph.getVertex(nextVertID)))
			_connectVertex(steinerTree, self._graph.getVertex(nextVertID), bridgeVertID, L1Dist)
		_removeOverlaps(steinerTree, L1Dist)
		self._graph = steinerTree


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
		return	totalWeight


	def plotGraph(self):
		numStartingVerts = len(self._startingVerts)
		for ((vertID, vertName), (vertX, vertY)) in zip(self._startingVerts, self._startingVertCoords):
			startingPlot = plot.scatter(x=vertX,y=vertY, c=(0,0,0), marker="s", linewidths=(8,))

		for (steinerPtID, steinerVert) in self._graph._verts.items():
			steinerPlot = plot.scatter(x=steinerVert.X, y=steinerVert.Y, c=(1,0,0), marker="o", linewidths=(4,))
			plot.annotate(s=steinerVert.name, xy=(steinerVert.X, steinerVert.Y), xytext=(-10,-10),
					textcoords="offset points", ha="right")

		for (edgeID, edgeObj) in self._graph.getEdges():
			srcVert = self._graph.getVertex(edgeObj.fromVert)
			destVert = self._graph.getVertex(edgeObj.toVert)
			plot.plot((srcVert.X, destVert.X), (srcVert.Y, destVert.Y), "k-")

		plot.legend((startingPlot, steinerPlot), ("Starting Vertices", "Steiner Points"))
		plot.show()


if __name__ == "__main__":
	A = Vertex(0, 15, "A")
	B = Vertex(5, 20, "B")
	C = Vertex(16, 24, "C")
	D = Vertex(20, 20, "D")
	E = Vertex(33, 25, "E")
	F = Vertex(23, 11, "F")
	G = Vertex(35, 7, "G")
	H = Vertex(25, 0, "H")
	I = Vertex(10, 3, "I")

	grid = SmarterPrimSolver([A,B,C,D,E,F,G,H,I])
	grid.addStationWeightFunc(lambda g,v: 1.2 * pow(g.getDegreeOfVertex(v), 1.5))
	grid.solve()
	grid.graphStats(printOutput=True)
	grid.plotGraph()
