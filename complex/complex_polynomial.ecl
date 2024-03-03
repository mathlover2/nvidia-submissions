:- module(complex_polynomial).
:- use_module(complex).
:- lib(ic).

:- export complex_polynomial/3.
:- export durand_kerner/4.

complex_polynomial(L, X, R) :-
        reverse(L, Lrev),
        complex_polynomial0(Lrev, X, R), !.

complex_polynomial0([], _, c(0,0)) :- !.
complex_polynomial0([H], _, H) :- !.
complex_polynomial0([H1, H2], X, R) :-
        complex_mul_add(H1, X, H2, R), !.
complex_polynomial0([H | T], X, R) :-
        complex_mul_add(H, X, Rn, R),
        complex_polynomial0(T, X, Rn), !.

:- inline(complex_mul_add/4).
complex_mul_add(c(Hr, Hi), c(Xr, Xi), c(Rr, Ri), c(RNr, RNi)) :-
        RNr $= Hr + Xr * Rr - Xi * Ri,
        RNi $= Hi + Xi * Rr + Xr * Ri, !.

%% The Durand-Kerner algorithm for calculating the complex roots of a
%% polynomial simultaneously. This is based on the following Python
%% code I wrote beforehand:
%% 
%% def polyfromcoeffs(coeffs): # Note: this is equivalent 
%%                             # to complex_polynomial above!
%%     n = len(coeffs)
%%     def func(x):
%%         return sum(x^i * c for (i,c) in zip(reversed(range(n)), coeffs))
%%     return func

%% def durand_kerner(init, coeffs, maxiter):
%%     f = polyfromcoeffs(coeffs)
%%     n = len(coeffs) - 1
%%     acc, newacc = [init^i for i in range(n)], []
%%     for _ in range(maxiter):
%%         for j in range(n):
%%             acc[j] = acc[j] - f(acc[j])/durand_kerner_helper(j, acc)
%%     return acc

%% def durand_kerner_helper(ind, arr):
%%     arr_ind = arr[ind]
%%     return prod(arr_ind - arr[j] for j in range(len(arr)) if j != ind)
%%

%% Since I have not set up a proper set of operators to clean up the
%% syntax of complex.ecl, the code I write here will be more lengthy.

%% The following is the initial calling function. It initializes the
%% array which will contain the roots, and then calls the subroutine
%% durand_kerner0. It corresponds roughly to everything prior to the
%% `for` loop in the second python function.

durand_kerner(InitVal, Coeffs, MaxIter, Res) :-
        length(Coeffs, Np), N is Np - 1,
        dim(Acc, [N]),
        powers(N, InitVal, Powers),
        (foreach(Pwr, Powers), for(I, 1, N), param(Acc)
         do subscript(Acc, [I], Pwr)),
        durand_kerner0(Acc, MaxIter, 1, Res, N, Coeffs).

%% This subroutine handles the array-manipulation portion of the
%% algorithm It corresponds more or less to the `for` loop in the
%% Python version, though the complex arithmetic is handled entirely
%% by the next subroutine.

durand_kerner0(Res, 0, _, Res, _, _).
durand_kerner0(Acc, I, J, Res, N, Coeffs) :-
        dim(NewAcc, [N]),
        (for(K, 1, N, 1), param(J, NewAcc, Coeffs, Acc)
        do subscript(NewAcc, [K], X),
           (K #\= J 
           -> subscript(Acc, [K], X)
           ;  durand_kerner1(J, Acc, X, Coeffs))),
        (J #= N 
        -> (Jnew is 1,
            Inew is I-1);
           (Jnew is J+1, Inew is I)),
        durand_kerner0(NewAcc, Inew, Jnew, Res, N, Coeffs).

%% This subroutine handles most of the the complex arithmetic used in
%% the algorithm. In the Python code, this is essentially the line

%% acc[j] = acc[j] - f(acc[j])/durand_kerner_helper(j, acc)

%% with the helper function in the Python code corresponding to
%% durand_kerner2, given after this.

durand_kerner1(J, Acc, X, Coeffs) :-
        subscript(Acc, [J], U),
        complex_polynomial(Coeffs, U, FU),
        complex_neg(FU, NFU),
        array_list(Acc, AccList),
        maplist(durand_kerner2(U), AccList, MulList),
        complex_prod(MulList, DenomUp),
        complex_recip(DenomUp, Denom),
        complex_mul(Denom, NFU, Dimin),
        complex_add(U, Dimin, X).

%% This final subroutine implements a small bit of logic. If the
%% second argument is *exactly* U, then we have run into our original
%% root, so we set the third argument to c(1,0) == 1 + 0*I to avoid
%% having it count in the multiplication. Otherwise, we let the third
%% argument be equal to the first argument subtracted from the second
%% argument. 

durand_kerner2(U, U, c(1,0)).
durand_kerner2(U, X, R) :-
        complex_neg(U, NU),
        complex_add(X, NU, NR),
        complex_neg(NR, R).
        

powers(N, InitVal, Powers) :-
        (fromto([1, c(1,0), [c(1,0)]], In, Out, [N, _, RPowers]), param(InitVal)
         do In = [I, Xp, T],
            complex_mul(InitVal, Xp, X),
            Ip is I + 1,
            Out = [Ip, X, [X | T]]),
        reverse(RPowers, Powers).