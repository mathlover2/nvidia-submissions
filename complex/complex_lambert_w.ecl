:- use_module(complex).
:- lib(ic).

:- inline(complex_lambert_w/2).

%% The complex Lambert $W$-function, defined as the multivalued
%% inverse of $Y = X*exp(X)$.
complex_lambert_w(X, Y) :-
        X c_eq Y * exp(Y).

:- inline(complex_r_lambert_w/2).

%% The complex $r$-Lambert $W$-function, defined as the multivalued
%% inverse of $Y = X*exp(X) + r*X$.
complex_r_lambert_w(R, X, Y) :-
        X c_eq Y * exp(Y) + R*Y.

