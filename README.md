# **FastMVBP**

This repository implements  our algorithm of listing top-k s-biplxes.

### Compile the code

```
$ cd Code
$ make clean
$ make
```

### Running example

Suppose we want to list the top-2 1-biplexes with at least 3 vertices in both left side and right side from "youtube.g".

```
$ ./FastMVBP  -f "../dataset/youtube.g" -s 1 -l 3 -r 3 -k 2
```

Then we get the following result:

```
start reading ../dataset/youtube.g
reading finish
Left vertexs:94238,Right vertexs:30087,Edges:293360
_lb_L:192,_ub_L3086,_lb_R:3,_ub_R64,graph_size:12958
_lb_L:3,_ub_L48,_lb_R:272,_ub_R471,graph_size:0
tot time = 9.87683 seconds
tot cnt = 2
|V|: 319 |L|: 316 |R|: 3
L: 70 192 46979 ...
R: 94518 94272 94426 
|V|: 319 |L|: 316 |R|: 3
L: 25352 70 192 ...
R: 94518 94272 94426 
```

For clarity, we have omitted some vertices in L.

The details of the usage of our code are as follows.

### Options

The most important command line options are:

```
-f <Path to the input file graph>
-g <three parameter needed as above>
-s <s-value>
-l <the lowerbound of the number of vertices in left side>
-r <the lowerbound of the number of vertices in right side>
-k <top-k, the number of s-biplex to be found>
```

The common format of options are as follows.

```
$ ./FastMVBP -f "<file_path>" -s <s-value> -l <left lowerbound> -r <right lowerbound> -k <k-number>
```

### Input format

For our code, we only support the graph format as follows:

- The first line is three non-negative integers to denote the number of total vertices, the number of  vertices in left side, and the number of edges respectively.
- Then, in the following n lines (n is the number of total vertices), the first integer is the index of each vertex, which starts from 0. And other integers for each line is all the neighbors of the  vertex corresponding to the line.
- Remarkably, for the graph G(L,R), the indexes for vertices of L in range from 0 to |L|-1, and the indexes for vertices of R in range from |L| to |L|+|R|-1.

The example of our format as follows:

```
9 5 10
0 5 6
1 5
2 6 8
3 6 7
4 7 8
5 0 1
6 0 2 3
7 3 4
8 2 4
```

### File transformation

By the way, we support the transforming edges list into the format our code needs, as follows:

```
$ ./FastMVBP -f "the graph file" -g n n_L E
```

For the three parameters of -g, n is the number of total vertices, n_L is the number of vertices in left side, and E is the number of edges. Remarkably, we also require that the indexes of vertices in left side is strictly smaller than right side.

We will generate a file named "example.g" as the transformation result.