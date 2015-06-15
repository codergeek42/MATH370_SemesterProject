2################################################################################
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
#
# These data structures and related functions/methods are based heavily on
# the csuf::Edge, csuf::Graph, and csuf::prim C++ code created by Dr. Kevin
# Wortman, also under the Apache 2.0 license, and available online at
#	http://mirror.thecodergeek.com/src/graph.h
#	http://mirror.thecodergeek.com/src/prim.h
################################################################################
# This implements the Vertex, Edge, and Graph types and their related methods.
################################################################################

class Vertex:
	__slots__ = ("_xc", "_yc", "_name")

	def __init__(self, x, y, name=""):
		self._xc = int(x)
		self._yc = int(y)
		self._name = name

	def __repr__(self):
		return "Vertex({x}, {y}, r'{name}')".format(x=self._xc, y=self._yc, name=self._name)


	def __str__(self):
		return "{name}({x}, {y})".format(x=self._xc, y=self._yc, name=self._name)


	@property
	def name(self):
		return self._name

	@name.setter
	def name(self, val):
		raise AttributeError("Vertex name is immutable.")


	@property
	def X(self):
		return self._xc

	@X.setter
	def X(self, val):
		raise AttributeError("X-coordinate is immutable.")


	@property
	def Y(self):
		return self._yc

	@Y.setter
	def Y(self, val):
		raise AttributeError("Y-coordinate is immutable.")


	def __eq__(self, other):
		return (self._xc == other._xc and self._yc == other._yc)

	def __ne__(self, other):
		return not (self == other)


## NB: Graph class will keep ID->vertex, ID->edge pairing maps
class Edge:
	__slots__ = ("_id", "_from", "_to", "_weight")

	def __init__(self, edgeID, vertexFromID, vertexToID, edgeWeight):
		self._from = vertexFromID
		self._to = vertexToID
		self._weight = edgeWeight
		self._id = edgeID


	def __repr__(self):
		return "Edge({id}, {src}, {dest}, {weight})".format(
				id=self._id, src=self._from, dest=self._to,
				weight=self._weight)

	@property
	def fromVert(self):
		return self._from

	@fromVert.setter
	def fromVert(self, val):
		raise AttributeError("Vertex is immutable.")


	@property
	def toVert(self):
		return self._to

	@toVert.setter
	def toVert(self, val):
		raise AttributeError("Vertex is immutable.")


	@property
	def weight(self):
		return self._weight

	@weight.setter
	def weight(self, val):
		raise AttributeError("Weight is immutable.")


	@property
	def edgeID(self):
		return self._id

	@edgeID.setter
	def edgeID(self, val):
		raise AttributeError("Edge ID is immutable.")


	def isIncident(self, otherVertexID):
		return (self._from == otherVertexID) or (self._to == otherVertexID)


	def oppositeVertex(self, vertexID):
		assert self.isIncident(vertexID), "Vertices not incident."

		if (self._from == vertexID):
			return self._to
		else:
			return self._from


class Graph:
	__slots__ = ("_adj", "_edges", "_verts")

	def __init__(self):
		self._adj = dict()
		self._edges = dict()
		self._verts = dict()

	def numVerts(self):
		return len(self._adj)


	def numEdges(self):
		return len(self._edges)


	def isVertex(self, vertexID):
		return (vertexID in self._verts)


	def isEdge(self, edgeID):
		return (edgeID in self._edges)


	def addVertex(self, vert):
		assert type(vert) is Vertex, "Not a vertex."

		vertID = max(self._verts)+1 if len(self._verts) > 0 else 0
		self._verts[vertID] = vert
		self._adj[vertID] = []
		return vertID


	def addEdgeWithVertsAndWeight(self, vertFromID, vertToID, edgeWeight):
		assert self.isVertex(vertFromID), "Not a vertex: {ID}".format(ID=vertFromID)
		assert self.isVertex(vertToID), "Not a vertex: {ID}".format(ID=vertToID)


		alreadyEdge = self.getEdgeBetween(vertFromID, vertToID)
		if not alreadyEdge:
			newID = max(self._edges)+1 if len(self._edges) > 0 else 0
			self._edges[newID] = Edge(newID, vertFromID, vertToID, edgeWeight)
			self._adj[vertToID].append(newID)
			self._adj[vertFromID].append(newID)
			return newID
		else:
			return alreadyEdge.edgeID



	def addEdge(self, otherEdge):
		return self.addEdgeWithVertsAndWeight(otherEdge.fromVert, otherEdge.toVert, otherEdge.weight)


	def getIncidentEdges(self, vertexID):
		assert self.isVertex(vertexID), "Not a vertex."

		return self._adj[vertexID]


	def getEdgeBetween(self, vertOneID, vertTwoID):
		assert self.isVertex(vertOneID), "Not a vertex."
		assert self.isVertex(vertTwoID), "Not a vertex."

		for edgeID in self.getIncidentEdges(vertOneID):
			if self._edges[edgeID].isIncident(vertTwoID):
				return self._edges[edgeID]
		return None

	def getNeighborhood(self, vertID):
		assert self.isVertex(vertID), "Not a vertex."

		neighborVerts = []
		for edgeID in self.getIncidentEdges(vertID):
			neighborVerts.append(self.getEdge(edgeID).oppositeVertex(vertID))

		return neighborVerts

	def getEdge(self, edgeID):
		assert self.isEdge(edgeID), "Not an edge."

		return self._edges[edgeID]

	def getVertex(self, vertID):
		assert self.isVertex(vertID), "Not a vertex."

		return self._verts[vertID]


	def getEdges(self):
		return list(self._edges.items())


	def getVerts(self):
		return list(self._verts.items())


	def removeEdge(self, edgeID):
		assert self.isEdge(edgeID), "Not an edge."


		for vertID in [self.getEdge(edgeID).fromVert, self.getEdge(edgeID).toVert]:
			self._adj[vertID].remove(edgeID)
		del self._edges[edgeID]


	def removeVertex(self, vertID):	
		assert self.isVertex(vertID), "Not a vertex."

		for edgeID in self._adj[vertID][:]:
			self.removeEdge(edgeID)
		del self._verts[vertID]


	def getDegreeOfVertex(self, vertID):
		return len(self._adj[vertID])


	def printGraph(self):
		print("Graph Vertices (total {n}):".format(n=len(self._verts)))
		for (vertID, vertObj) in self._verts.items():
			print("\t({ID}) {vertInfo!s}".format(ID=vertID, vertInfo=vertObj))

		print("Graph Edges (total {m})".format(m=len(self._edges)))
		totalWeight = 0
		for (edgeID, edgeObj) in self._edges.items():
			print("\t({ID}) {src!s} -- {dest!s}  [{weight:.2f}]".format(
				ID=edgeID, src=self._verts[edgeObj.fromVert], 
				dest=self._verts[edgeObj.toVert], weight=edgeObj.weight))
			totalWeight += edgeObj.weight

		print("Total weight: {weight:.2f}".format(weight=totalWeight))


	def getVertexID(self, x, y):
		for (vertID, vertObj) in self._verts.items():
			if vertObj == Vertex(x, y):
				return vertID
		return -1


	def getSubgraph(self, edgeList):
		subgraph = Graph()
		subgraph._verts = self._verts
		subgraph._adj = dict([(vertID, []) for vertID in subgraph._verts.keys()]) 
		for edgeID in sorted(edgeList):
			subgraph.addEdge(self._edges[edgeID])
		return subgraph


	def edgeStr(self, edgeID):
		assert self.isEdge(edgeID), "Not an edge."

		edge = self._edges[edgeID]
		return "{src!s} --> {dest!s} [{weight:.5f}]".format(src=self._verts[edge.fromVert],
				dest=self._verts[edge.toVert], weight=edge.weight)



def L1Dist(vert1, vert2):
	return abs(vert2.X-vert1.X) + abs(vert2.Y-vert1.Y)



def pythagDist(vert1, vert2):
	from math import sqrt

	return sqrt((vert2.X - vert1.X)**2 + (vert2.Y - vert1.Y)**2)


def makeGraphComplete(theGraph, weightFunc):
	from itertools import combinations

	allVerts = theGraph.getVerts()
	allVertPairs = combinations(allVerts, 2)
	for (vert1, vert2) in allVertPairs:
		(vertFromID, vertFrom) = vert1
		(vertToID, vertTo) = vert2
		theGraph.addEdgeWithVertsAndWeight(vertFromID, vertToID, weightFunc(vertFrom, vertTo))


