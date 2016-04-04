from __future__ import generators

def minKey(key, mstSet):
    min1 = float("inf")
    for v in range(n):
        if mstSet[v] == False and key[v] < min1:
            min1 = key[v]
            minInd = v
    return minInd

def printMST(par, n, mat):
    print("Printing MST")
    print("Edge Weight")
    for i in range(n):
        if par[i] != -1:
            print(par[i], i, mat[i][par[i]])
        
def mstPrim(mat, start):
    print("\n\n\n")
    par = [0  for i in range(n)]
    key = [float("inf") for i in range(n)]
    mstSet = [False for i in range(n)]
    key[start] = 0
    par[start] = -1

    for c in range(1,n):
        u = minKey(key, mstSet)
        mstSet[u] = True
        for v in range(n):
            if mat[u][v] and mstSet[v] == False and mat[u][v] < key[v]:
                par[v] = u
                key[v] = mat[u][v]
    printMST(par, n, mat)

if 'True' not in globals():
    globals()['True'] = not None
    globals()['False'] = not True

class unionFind:
    def __init__(self):
        self.weights = {}
        self.parents = {}

    def __getitem__(self, object):
        if object not in self.parents:
            self.parents[object] = object
            self.weights[object] = 1
            return object
        
        path = [object]
        root = self.parents[object]
        while root != path[-1]:
            path.append(root)
            root = self.parents[root]
        
        for ancestor in path:
            self.parents[ancestor] = root
        return root

    def union(self, *objects):
        roots = [self[x] for x in objects]
        heaviest = max([(self.weights[r],r) for r in roots])[1]
        for r in roots:
            if r != heaviest:
                self.weights[heaviest] += self.weights[r]
                self.parents[r] = heaviest

def matching(G, initialMatching = {}):
    matching = {}
    for x in initialMatching:
        matching[x] = initialMatching[x]

    for v in G:
        if v not in matching:
            for w in G[v]:
                if w not in matching:
                    matching[v] = w
                    matching[w] = v
                    break

    def augment():
        leader = unionFind()
        S = {}
        T = {}
        unexplored = []
        base = {}
        
        def blossom(v,w,a):
            def findSide(v,w):
                path = [leader[v]]
                b = (v,w)
                while path[-1] != a:
                    tnode = S[path[-1]]
                    path.append(tnode)
                    base[tnode] = b
                    unexplored.append(tnode)
                    path.append(leader[T[tnode]])
                return path
            
            a = leader[a]
            path1,path2 = findSide(v,w), findSide(w,v)
            leader.union(*path1)
            leader.union(*path2)
            S[leader[a]] = S[a]

        topless = object()
        def alternatingPath(start, goal = topless):
            path = []
            while 1:
                while start in T:
                    v, w = base[start]
                    vs = alternatingPath(v, start)
                    vs.reverse()
                    path += vs
                    start = w
                path.append(start)
                if start not in matching:
                    return path
                tnode = matching[start]
                path.append(tnode)
                if tnode == goal:
                    return path
                start = T[tnode]
                
        def pairs(L):
            i = 0
            while i < len(L) - 1:
                yield L[i],L[i+1]
                i += 2
            
        def alternate(v):
            path = alternatingPath(v)
            path.reverse()
            for x,y in pairs(path):
                matching[x] = y
                matching[y] = x

        def addMatch(v, w):
            alternate(v)
            alternate(w)
            matching[v] = w
            matching[w] = v
            
        def ss(v,w):
            if leader[v] == leader[w]:
                return False       
            path1, head1 = {}, v
            path2, head2 = {}, w
    
            def step(path, head):
                head = leader[head]
                parent = leader[S[head]]
                if parent == head:
                    return head
                path[head] = parent
                path[parent] = leader[T[parent]]
                return path[parent]
                
            while 1:
                head1 = step(path1, head1)
                head2 = step(path2, head2)
                
                if head1 == head2:
                    blossom(v, w, head1)
                    return False
                
                if leader[S[head1]] == head1 and leader[S[head2]] == head2:
                    addMatch(v, w)
                    return True
                
                if head1 in path2:
                    blossom(v, w, head1)
                    return False
                
                if head2 in path1:
                    blossom(v, w, head2)
                    return False    

        for v in G:
            if v not in matching:
                S[v] = v
                unexplored.append(v)

        current = 0
        while current < len(unexplored):
            v = unexplored[current]
            current += 1

            for w in G[v]:
                if leader[w] in S:
                    if ss(v,w):
                        return True

                elif w not in T:
                    T[w] = v
                    u = matching[w]
                    if leader[u] not in S:
                        S[u] = w
                        unexplored.append(u)
                        
        return False
                        
    while augment():
        pass

    return matching

def edmonds_karp(C, source, sink):
    n = len(C)
    F = [[0] * n for i in range(n)]

    while True:
        path = bfs(C, F, source, sink)
        if not path:
            break
        flow = min(C[u][v] - F[u][v] for u,v in path)
        for u, v in path:
            F[u][v] += flow
            F[v][u] -= flow
        print(F)
    return sum(F[source][i] for i in range(n))

def bfs(C, F, source, sink):
    queue = [source]                 
    paths = {source: []}
    while queue:
        u = queue.pop(0)
        for v in range(len(C)):
            if C[u][v] - F[u][v] > 0 and v not in paths:
                paths[v] = paths[u] + [(u, v)]
                if v == sink:
                    return paths[v]
                queue.append(v)
    return None

n = int(input("Enter the number of vertices : "))

mat = [[float("inf") for i in range(n)] for i in range(n)]

numEdg = int(input("Enter number of edges : "))

print("Enter edges in u, v, w format where the edge is between vertices u, v with weight w")

for i in range(1, numEdg + 1):
    u, v, w = map(int, input("Enter the edge details : ").split())
    mat[u][v] = w
    mat[v][u] = w

print("Prim's algorithm")
start = int(input("Enter the starting vertex : "))
mstPrim(mat, start)

print("\n\n\nEdmond's algorithm")
G = {}
for z in range(n):
    app = []
    for p in range(n):
        if mat[z][p] > 0 and mat[z][p] < float("inf"):
            app.append(p)
    G[z] = app

print("Edges")    
print(matching(G))
