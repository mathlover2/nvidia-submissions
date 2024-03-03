:- module(complex).
:- lib(ic).
:- import nth1/3 from listut.
:- export struct(c(re, im)).

:- export re/2, im/2, complex_add/3, complex_mul/3, complex_neg/2.
:- export complex_conjugate/2, complex_recip/2, complex_exp/2,
   complex_abs/2.
:- export complex_sqrt/2, complex_sinh/2, complex_sin/2.
:- export complex_cosh/2, complex_cos/2, complex_log/2, complex_log/3.
:- export complex_atanh/2, complex_atan/2, complex_tan/2.
:- export complex_tanh/2, complex_asinh/2, complex_acosh/2,
   complex_asin/2.
:- export complex_acos/2, atan2_raw/3.
:- export complex_eq/2, complex_sum/2, complex_prod/2, complex_pow/3.
:- export abs_leq/2, abs_geq/2, c_locate/2, c_locate/3, sqr/2.

:- local  op(700, xfx, $:=).
:- export op(700, xfy, [c_eq, c_is, abs_leq, abs_geq]).
:- export op(200, xfx, c).
:- export op(210, fy, /).

:- export c_eq/2, c_is/2.

:- inline(re/2).
:- inline(im/2).

re(c(R, _), R).
im(c(_, I), I).


complex_eq(c(R1, R2), c(I1, I2)) :-
    R1 $:= I1,
    R2 $:= I2.

r_to_c(R, c(R, 0)).


complex_add(c(R1, I1), c(R2, I2), c(R3, I3)) :-
    R3 $:= R1 + R2,
    I3 $:= I1 + I2.

complex_mul(c(R1, I1), c(R2, I2), c(R3, I3)) :-
    R3 $:= R1 * R2 - I1 * I2,
    I3 $:= R2 * I1 + R1 * I2, !.

complex_neg(c(R1, I1), c(R2, I2)) :-
    R2 $:= -R1,
    I2 $:= -I1.

complex_conjugate(c(R1, I1), c(R1, I2)) :-
    I2 $:= -I1.

complex_recip(c(R1, I1), c(R2, I2)) :-
    Norm $:= sqr(R1) + sqr(I1),
    R2 $:= R1 / Norm,
    I2 $:= -I1 / Norm.


complex_exp(c(R, I), Z2) :-
    ExpR $:= exp(R),
    RotR $:= cos(I),
    RotI $:= sin(I),
    complex_mul(c(ExpR, 0), c(RotR, RotI), Z2).


complex_abs(c(R, I), Abs) :-
    Abs $:= sqrt(sqr(R) + sqr(I)).

complex_sqrt(c(R, 0), Z2) :-
    R $>= 0, SqrtR $:= sqrt(R), Z2 = c(SqrtR, 0), !.

complex_sqrt(c(R, 0), Z2) :-
    R $< 0, SqrtR $:= sqrt(-R), Z2 = c(0, SqrtR), !.

complex_sqrt(c(R, I), Z2) :-
    complex_abs(c(R, I), Abs),
    SqrtAbs $:= sqrt(Abs),
    atan2_raw(I, R, Td), T $:= Td/2,
    complex_exp(c(0, T), ExpT),
    complex_mul(c(SqrtAbs, 0), ExpT, Z2).


complex_sinh(Z1, Z2) :-
    OneHalf = c(0.5, 0),
    complex_neg(Z1, NegZ1),
    complex_exp(NegZ1, ExpNegZ1),
    complex_neg(ExpNegZ1, NegExpNegZ1),
    complex_exp(Z1, ExpZ1),
    complex_add(ExpZ1, NegExpNegZ1, Twice),
    complex_mul(Twice, OneHalf, Z2).


complex_sin(Z1, Z2) :-
    I = c(0, 1),
    complex_mul(Z1, I, Sl1),
    complex_sinh(Sl1, Sl2),
    complex_mul(Sl2, I, Sl3),
    complex_neg(Sl3, Z2).


complex_cosh(Z1, Z2) :-
    OneHalf = c(0.5, 0),
    complex_neg(Z1, NegZ1),
    complex_exp(NegZ1, ExpNegZ1),
    complex_exp(Z1, ExpZ1),
    complex_add(ExpZ1, ExpNegZ1, Twice),
    complex_mul(Twice, OneHalf, Z2).


complex_cos(Z1, Z2) :-
    I = c(0, 1),
    complex_mul(Z1, I, Sl),
    complex_cosh(Sl, Z2).


complex_log(Z1, Z2) :-
    complex_log(Z1, Z2, 0).


%% Principal branch of the log function.

complex_log(c(R, I), Z2, K) :-
        I $= 0,
        LogR $:= ln(abs(R)),
        (R $> 0 -> T1 $:= 0; T1 $:= pi),
        T $:= T1 + 2 * pi * K,
        c(LogR, T) = Z2.

complex_log(c(R, I), Z2, K) :-
        R $= 0,
        LogR $:= ln(abs(I)),
        (I $> 0 -> T1 $:= pi/2 ; T1 $:= -pi/2),
        T $:= T1 + 2 * pi * K,
        c(LogR, T) = Z2.
        
complex_log(Z1, Z2, K) :-
    c(R, I) = Z1,
    complex_abs(Z1, Rr),
    atan2_raw(I, R, T1),
    T $:= T1 + 2 * pi * K,
    LogR $:= ln(Rr),
    c(LogR, T) = Z2.


complex_atanh(Z1, Z2) :-
    complex_atanh(Z1, Z2, 0).


complex_atanh(Z1, Z2, K) :-
    complex_add(c(1, 0), Z1, Sl1),
    complex_add(c(-1, 0), Z1, Sl2),
    complex_neg(Sl2, Sl3),
    complex_recip(Sl3, Sl4),
    complex_mul(Sl1, Sl4, Res1),
    complex_log(Res1, Res2, K),
    complex_mul(c(0.5, 0), Res2, Z2).


complex_atan(Z1, Z2) :-
    complex_atan(Z1, Z2, 0).


complex_atan(Z1, Z2, K) :-
    complex_mul(Z1, c(0, 1), Sl),
    complex_atanh(Sl, Sl2, K),
    complex_mul(c(0, -1), Sl2, Z2).


complex_tan(Z1, Z2) :-
    complex_sin(Z1, SinZ1),
    complex_cos(Z1, CosZ1),
    complex_recip(CosZ1, CosZ1Recip),
    complex_mul(CosZ1Recip, SinZ1, Z2).


complex_tanh(Z1, Z2) :-
    complex_sinh(Z1, SinhZ1),
    complex_cosh(Z1, CoshZ1),
    complex_recip(CoshZ1, CoshZ1Recip),
    complex_mul(CoshZ1Recip, SinhZ1, Z2).


complex_asinh(Z1, Z2) :-
    complex_sqr(Z1, Z1Square),
    One c_is 1,
    complex_add(Z1Square, One, Loc1),
    complex_sqrt(Loc1, Loc2),
    complex_add(Loc2, Z1, Loc3),
    complex_log(Loc3, Z2).


complex_acosh(Z1, Z2) :-
    complex_mul(Z1, Z1, Z1Square),
    complex_add(Z1Square, c(-1, 0), Loc1),
    complex_sqrt(Loc1, Loc2),
    complex_add(Loc2, Z1, Loc3),
    complex_log(Loc3, Z2).


complex_asin(Z1, Z2) :-
    complex_mul(c(0, 1), Z1, Loc),
    complex_asinh(Loc, Loc1),
    complex_mul(c(0, -1), Loc1, Z2).


complex_acos(Z1, Z2) :-
    HalfPi $:= pi / 2,
    complex_asin(Z1, Loc),
    complex_neg(Loc, Loc1),
    complex_add(c(HalfPi, 0), Loc1, Z2).
        
atan2_raw(Y, X, Res) :-
        (X $> 0 and Res $= atan(Y/X))
        or (X $< 0 and Y $>= 0 and Res $= atan(Y/X) + pi)
        or (X $< 0 and Y $< 0 and Res $= atan(Y/X) - pi)
        or (X $= 0 and Y $> 0 and Res $= pi/2)
        or (X $= 0 and Y $< 0 and Res $= -pi/2)
        or (X $= 0 and Y $= 0 and Res $= 0).


%% More complicated functions built on the basic ones defined above.

complex_sum([], c(0, 0)).
complex_sum([H], H).
complex_sum(Ss, c(Rt, It)) :-
        (foreach(R, Rs), foreach(I, Is), foreach(S, Ss)
        do S = c(R, I)),
        Rt $:= sum(Rs), It $:= sum(Is).

complex_prod([], 1 c 0).
complex_prod(L, Res) :-
        complex_prod0(L, c(1, 0), c(R_expr, I_expr)),
        R $:= eval(R_expr),
        I $:= eval(I_expr),
        Res = c(R, I).

complex_prod0([], Res, Res).
complex_prod0([c(R0, I0) | T], c(Ra, Ia), Res) :-
        Rn = c(R0 * Ra - I0 * Ia, Ra * I0 + R0 * Ia),
        complex_prod0(T, Rn, Res).


%% Complex powers with integer exponents.

complex_pow(_, 0, c(1, 0)).
complex_pow(X, 1, X) :- !.
complex_pow(X, 2, Y) :-
        complex_sqr(X, Y), !.
complex_pow(X, 3, Y) :-
        complex_cube(X, Y), !.
complex_pow(X, 4, Y) :-
       complex_fourth(X, Y), !.
complex_pow(X, 5, Y) :-
       complex_fifth(X, Y), !.
complex_pow(X, 6, Y) :-
       complex_sixth(X, Y), !.
complex_pow(X, N, Y) :-
        integer(N), N < 0,
        Xr c_is / X,
        NegN is -N,
        complex_pow(Xr, NegN, Y).

%% This version of the function calculates the real part and the
%% imaginary part of the result via the binomial theorem. It uses
%% two features to optimize the resulting system of constraints:
%%
%% 1. It calculates the real and imaginary parts independently. These
%% auxiliary functions can be found later.

%% 2. It explicitly reuses intermediate values via lists of powers of
%% the real and imaginary parts of the first argument. This is done
%% using the power_list/3 predicate, given later as well.

complex_pow(X1 c X2, N, Y1 c Y2) :-
        integer(N), N > 3, 
        pow_list(X1, N, X1_powlist),
        pow_list(X2, N, X2_powlist),
        complex_pow_real(N, X1_powlist, X2_powlist, Y1),
        complex_pow_imag(N, X1_powlist, X2_powlist, Y2).

%% Complex powers with complex exponents.
complex_pow(X, Z, Y) :- 
        Z = _ c _,
        complex_log(X, LogX),
        complex_mul(LogX, Z, ZLogX),
        complex_exp(ZLogX, Y).

$:=(X, Y) :-
        float(X), float(Y), X =:= Y, !.
$:=(X, Y) :-
        var(X), var(Y), X = Y, !.
$:=(X, Y) :-
        var(X), ground(Y),
        X is Y, !.
$:=(X, Y) :-
        eval(X) $= eval(Y).

c_is(c(X1, X2), c(X3, X4)) :-
        X1 $:= X3,
        X2 $:= X4,
        !.

c_is(c(X, 0), X) :- number(X), !.
c_is(c(X, 0), X) :- var(X), !.

c_is(X, Ye + Ze) :-
        Y c_is Ye,
        Z c_is Ze,
        complex_add(Y, Z, X).

c_is(X, - Ye) :-
        Y c_is Ye,
        complex_neg(Y, X).

c_is(X, Ye - Ze) :-
        Y c_is Ye,
        Zneg c_is -Ze,
        complex_add(Y, Zneg, X).

c_is(X, Ye * Ze) :-
        Y c_is Ye,
        Z c_is Ze,
        complex_mul(Y, Z, X).

c_is(X, / Ye) :-
        Y c_is Ye,
        complex_recip(Y, X).

c_is(X, Ye / Ze) :-
        Y c_is Ye,
        Z c_is Ze,
        complex_recip(Z, Zr),
        complex_mul(Y, Zr, X).

c_is(X, Ye ^ Ze) :-
        (Ze = _ c _; real(Ze)),
        Y c_is Ye,
        Z c_is Ze,
        complex_pow(Y, Z, X).

c_is(X, Ye ^ Ne) :-
        integer(Ne),
        Y c_is Ye,
        N is Ne,
        complex_pow(Y, N, X).


c_is(X, \ Ye) :-
        Y c_is Ye,
        complex_conjugate(Y, X).

c_is(c(X, 0), abs(Ye)) :-
        Y c_is Ye,
        complex_abs(Y, X).

c_is(X, exp(Ye)) :-
        Y c_is Ye,
        complex_exp(Y, X).

c_is(X, sqr(Ye)) :-
        Y c_is Ye,
        complex_sqr(Y, X).

c_is(X, cube(Ye)) :-
        Y c_is Ye,
        complex_cube(Y, X).

c_is(X, sqrt(Ye)) :-
        Y c_is Ye,
        complex_sqrt(Y, X).

c_is(X, sin(Ye)) :-
        Y c_is Ye,
        complex_sin(Y, X).

c_is(X, sinh(Ye)) :-
        Y c_is Ye,
        complex_sinh(Y, X).

c_is(X, cos(Ye)) :-
        Y c_is Ye,
        complex_cos(Y, X).

c_is(X, cosh(Ye)) :-
        Y c_is Ye,
        complex_cosh(Y, X).

c_is(X, log(Ye)) :-
        Y c_is Ye,
        complex_log(Y, X).

c_is(X, ln(Ye)) :-
        Y c_is Ye,
        complex_log(Y, X).

c_is(X, tan(Ye)) :-
        Y c_is Ye,
        complex_tan(Y, X).

c_is(X, tanh(Ye)) :-
        Y c_is Ye,
        complex_tanh(Y, X).

c_is(X, atan(Ye)) :-
        Y c_is Ye,
        complex_atan(Y, X).

c_is(X, atanh(Ye)) :-
        Y c_is Ye,
        complex_atanh(Y, X).

c_is(X, asin(Ye)) :-
        Y c_is Ye,
        complex_asin(Y, X).

c_is(X, asinh(Ye)) :-
        Y c_is Ye,
        complex_asinh(Y, X).

c_is(X, acos(Ye)) :-
        Y c_is Ye,
        complex_acos(Y, X).

c_is(X, acosh(Ye)) :-
        Y c_is Ye,
        complex_acosh(Y, X).

c_is(X, asin(Ye)) :-
        Y c_is Ye,
        complex_asin(Y, X).

c_is(X, sum(Le)) :-
        (foreach(Le0, Le), foreach(S0, S) do S0 c_is Le0),
        complex_sum(S, X).

c_is(X, prod(Le)) :-
        (foreach(Le0, Le), foreach(S0, S) do S0 c_is Le0),
        complex_prod(S, X).

c_eq(Xe, Ye) :-
        X c_is Xe,
        Y c_is Ye,
        X c_is Y.

abs_leq(Xe, B) :-
        X1 c X2 c_is Xe,
        sqr(B) $>= sqr(X1) + sqr(X2), !.

abs_leq([], _) :- !.
abs_leq([X0e | XeT], B) :-
        abs_leq(X0e, B), !,
        abs_leq(XeT, B).

abs_geq(Xe, B) :-
        X1 c X2 c_is Xe,
        sqr(B) $=< sqr(X1) + sqr(X2), !.

abs_geq([], _) :- !.
abs_geq([X0e | XeT], B) :-
        abs_geq(X0e, B), !,
        abs_geq(XeT, B).

c_locate(CVars, Prec) :-
        extract_vars(CVars, Vars),
        locate(Vars, Prec).

c_locate(CVars, Prec, LinLog) :-
        extract_vars(CVars, Vars),
        locate(Vars, Prec, LinLog).


extract_vars([], []).
extract_vars([V1 c V2 | CVt], [V1, V2 | Vt]) :-
        extract_vars(CVt, Vt).
        

%% Helper functions.

complex_sqr(X1 c X2, Y1 c Y2) :-
        Y1 $:= sqr(X1) - sqr(X2),
        Y2 $:= 2*X1*X2.

complex_cube(X1 c X2, Y1 c Y2) :-
        %% Reusing these values allows us to save on the number of
        %% delayed goals a little.
        X1_2 $:= sqr(X1),
        X2_2 $:= sqr(X2),
        
        Y1 $:= (X1_2 - 3*X2_2)*X1,
        Y2 $:= (3*X1_2 - X2_2)*X2.

complex_fourth(X1 c X2, Y1 c Y2) :-
	A2 $:= sqr(X1),
	B2 $:= sqr(X2),
        AB $:= X1*X2,
	Y1 $:= sqr(A2) - 6*sqr(AB) + sqr(B2),
	Y2 $:= 4*(A2 - B2)*AB.

complex_fifth(X1 c X2, Y1 c Y2) :-
	A2 $:= sqr(X1),
	B2 $:= sqr(X2),
        AB $:= X1*X2,
	Y1 $:= (5*sqr(B2) - 10*sqr(AB) + sqr(A2))*X1,
	Y2 $:= (5*sqr(A2) - 10*sqr(AB) + sqr(B2))*X2.

complex_sixth(X, Y) :-
	Y c_eq sqr(X^3).

complex_pow_real(N, X1_powlist, X2_powlist, Y1) :-
        (for(I, 0, N, 2), fromto(Expr, In, Out, _), param(N, X1_powlist, X2_powlist)
        do complex_pow_construct_term(N, I, X1_powlist, X2_powlist, Term, true),
           ((I =:= N; I =:= N - 1) -> In = Term ; In = Term + Out)),
        Y1 $= eval(Expr).

complex_pow_imag(N, X1_powlist, X2_powlist, Y2) :-
        (for(I, 1, N, 2), fromto(Expr, In, Out, _), param(N, X1_powlist, X2_powlist)
        do complex_pow_construct_term(N, I, X1_powlist, X2_powlist, Term, false),
           ((I =:= N; I =:= N - 1) -> In = Term ; In = Term + Out)),
        Y2 $= eval(Expr).

complex_pow_construct_term(N, I, X1_powlist, X2_powlist, Term, IsReal) :-
        binom(N, I, C1),
        ((IsReal -> (I =:= 0 -> nth1(N, X1_powlist, Term));
                    (I =:= 1 -> (Ni is N - 1, 
                                 nth1(Ni, X1_powlist, A),
                                 nth1(1, X2_powlist, B),
                                 Term = N*A*B)))
        ; (I =:= N - 1, 1 is (-1)^(I // 2)) 
        -> (nth1(1, X1_powlist, A), nth1(I, X2_powlist, B), Term = N*A*B)
        ; (I =:= N - 1, -1 is (-1)^(I // 2)) 
        -> (nth1(1, X1_powlist, A), nth1(I, X2_powlist, B), Term = -N*A*B)
        ; (I =:= N, 1 is (-1)^(I // 2)) -> nth1(N, X2_powlist, Term)
        ; (I =:= N, -1 is (-1)^(I // 2)) -> (nth1(N, X2_powlist,
                                                  B), Term = -B)
        ; (C is C1 * (-1)^(I // 2),
           Idx is N - I,
           nth1(Idx, X1_powlist, A),
           nth1(I, X2_powlist, B),
           Term = C * A * B)), !.

pow_list(X, N, List) :-
        (for(I, 1, N), foreach(L, List), param(X)
        do (I =:= 1 -> L = X; I =:= 2 -> L $:= sqr(X); L $:= X^I)).

binom(N, K, R) :-
        (for(I, 1, K, 1),
         fromto([1, 1], [In1, In2],
                [Out1, Out2], [Num, Den]), param(N)
        do  Out1 is In1 * (N - I + 1),
            Out2 is In2 * I),
        R is Num // Den.

sqr(X, Y) :-
        Y is X*X.