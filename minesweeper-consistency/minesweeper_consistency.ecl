:- lib(ic).
:- import occurrences/3 from ic_global.
:- import flatten/2 from lists.

% Constraints

minesweeper_consistency(Grid) :-
        (foreachelem(X, Grid, Idx), param(Grid)
        do neighbor_constraint(X, Grid, Idx)).

neighbor_constraint(X, Grid, Idx) :-
        neighbors(Grid, Idx, Neighbors, NL),
        NLm is NL - 1,
        X #:: [-1..NLm],
        BombCount #:: [0..NLm],
        occurrences(-1, Neighbors, BombCount),
        (X #= BombCount) or (X #= -1).

% Gives the neighbors of the entry given by Idx in Grid.

neighbors(Grid, Idx, Neighbors, NL) :-
        Idx = [Y, X],
        dim(Grid, [N, M]),
        (X == 1 -> NeighborXIndices = [1,2];
         X == M -> (Mm is M - 1, NeighborXIndices = [Mm, M]);
         (Xm is X-1, Xp is X+1, NeighborXIndices = [Xm, Xp])),
        (Y == 1 -> NeighborYIndices = [1,2];
         Y == N -> (Nm is N - 1, NeighborYIndices = [Nm, N]);
         (Ym is Y-1, Yp is Y+1, NeighborYIndices = [Ym, Yp])),
        [X1, X2] = NeighborXIndices,
        [Y1, Y2] = NeighborYIndices,
        subscript(Grid, [Y1..Y2, X1..X2], NeighborsArray),
        eval_to_array(concat(NeighborsArray), Neighbors),
        dim(Neighbors, [NL]).

% Setup problem from database.

setup_problem(Number, Matrix) :-
        problem(Number, List, Dims),
        flatten(List, Lc),
        list_matrix(Lc, Dims, Matrix),
        minesweeper_consistency(Matrix).

%% Search predicate

minesweeper_search(Matrix) :-
        search(Matrix, 0, input_order, choice, complete, []).

choice(X) :- X #>= 0.
choice(X) :- X #= -1.

% Display a matrix in a nicer form.

show_matrix(Matrix) :-
        array_list(Matrix, List1),
        maplist(array_list, List1, OutputList),
        maplist(writeln, OutputList).
        

% Helper functions.

list_matrix(List, [M, N], Matrix) :-
        dim(Matrix, [M, N]),
        eval_to_list(concat(Matrix), List).

% Sample problem database

problem(1, [[0,0,0,0,0,1,_,_],
            [0,0,0,0,0,1,2,_],
            [0,0,0,0,0,0,1,_],
            [0,0,0,0,0,0,2,_],
            [0,0,0,0,0,0,1,_],
            [2,3,2,1,1,2,3,_],
            [_,_,_,_,_,_,_,_]], [7, 8]).

problem(2, [[_,_,_,_,_,_,_],
            [_,2,1,1,2,1,1],
            [_,2,0,0,0,0,0],
            [_,2,1,0,0,0,0],
            [_,_,1,0,0,0,0],
            [_,_,2,0,0,1,1],
            [_,_,1,0,0,1,_]], [7, 7]).

problem(3, [[F,D,2,1,2,1],
            [A,A,3,-1,4,B],
            [2,2,3,-1,5,B],
            [0,0,1,1,4,B],
            [0,1,1,1,2,B],
            [0,1,C,E,E,E]], [6, 6]).

problem(4, [[1,1,1,0],
            [1,-1,2,1],
            [1,2,_,_],
            [0,1,_,_]], [4, 4]).

problem(5, [[0,0,0,1,1],
            [0,0,0,2,_],
            [0,0,0,2,_],
            [1,1,0,1,1],
            [_,1,0,0,0]], [5, 5]).

problem(6, [[0,_,1],
            [2,_,5],
            [_,_,_]], [3, 3]).

problem(7, [[_,_,_,_,_,_,_,_,_,_],
            [_,_,_,_,_,_,_,_,2,_],
            [1,2,_,1,1,1,2,1,1,_],
            [0,1,1,1,0,0,0,0,1,_],
            [0,0,0,0,0,0,0,1,2,_],
            [1,1,0,0,0,0,0,1,_,_],
            [_,2,0,1,2,2,1,1,_,_],
            [_,3,1,1,_,_,_,_,_,_],
            [_,_,_,_,_,_,_,_,_,_]], [9, 10]).

problem(8, [[_,_,_,_,_,_],
            [_,2,2,2,2,_],
            [_,2,0,0,2,_],
            [_,2,0,0,2,_],
            [_,2,2,2,2,_],
            [_,_,_,_,_,_]], [6, 6]).

problem(9, [[2,3,-1,2,2,-1,2,1],
            [-1,-1,5,_,_,4,-1,2],
            [_,_,-1,-1,-1,_,4,-1],
            [-1,6,_,6,-1,-1,_,2],
            [2,-1,-1,_,5,5,_,2],
            [1,3,4,_,-1,-1,4,-1],
            [0,1,-1,4,_,_,-1,3],
            [0,1,2,-1,2,3,-1,2]], [8, 8]).

problem(10, [[_,2,0,0,0,1,_],
             [_,2,0,0,1,2,_],
             [_,2,1,0,1,_,_],
             [_,_,2,0,1,_,_],
             [_,_,3,1,1,_,_],
             [_,_,_,_,_,_,_]], [6, 7]).

problem(11, [[_,_,_,_,_,_],
             [_,_,2,1,2,_],
             [1,1,1,0,1,_],
             [0,0,0,0,2,_],
             [0,1,1,1,2,_],
             [0,1,_,_,_,_],
             [1,2,_,_,_,_],
             [_,_,_,_,_,_]], [8, 6]).

problem(12, [[2,-1,_],
             [2,-1,_], 
             [1,3,_],
             [0,2,_],
             [0,2,_],
             [0,1,_],
             [0,2,_],
             [0,1,_]], [8, 3]).