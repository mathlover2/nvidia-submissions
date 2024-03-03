:- module(complex_noic).
:- local struct(c(re, im)).

:- export re/2, im/2, complex_add/3, complex_mul/3, complex_neg/2.
:- export complex_conjugate/2, complex_recip/2, complex_exp/2,
   complex_abs/2.
:- export complex_sqrt/2, complex_sinh/2, complex_sin/2.
:- export complex_cosh/2, complex_cos/2, complex_log/2, complex_log/3.
:- export complex_atanh/2, complex_atan/2, complex_tan/2.
:- export complex_tanh/2, complex_asinh/2, complex_acosh/2,
   complex_asin/2.
:- export complex_acos/2, atan2_raw/3.
:- export complex_eq/2, complex_sum/2, complex_prod/2.

:- export op(700, xfy, c_is).
:- export op(200, xfx, c).
:- export op(210, fy, /).
:- export c_is/2.

:- inline(re/2).
:- inline(im/2).
re(c(R, _), R).
im(c(_, I), I).

:- inline(complex_eq/2).
complex_eq(c(R1, R2), c(I1, I2)) :-
        R1 is I1,
        R2 is I2.

:- inline(r_to_c/2).
r_to_c(R, c(R, 0)).

:- inline(complex_add/3).
complex_add(c(R1, I1), c(R2, I2), c(R3, I3)) :-
    R3 is R1 + R2,
    I3 is I1 + I2.

:- inline(complex_mul/3).
complex_mul(c(R1, I1), c(R2, I2), c(R3, I3)) :-
    R3 is R1 * R2 - I1 * I2,
    I3 is R2 * I1 + R1 * I2.

:- inline(complex_neg/2).
complex_neg(c(R1, I1), c(R2, I2)) :-
    R2 is -R1,
    I2 is -I1.

:- inline(complex_conjugate/2).
complex_conjugate(c(R1, I1), c(R1, I2)) :-
    I2 is -I1.

:- inline(complex_recip/2).
complex_recip(c(R1, I1), c(R2, I2)) :-
    Norm is R1^2 + I1^2,
    R2 is R1 / Norm,
    I2 is -I1 / Norm.

:- inline(complex_exp/2).
complex_exp(c(R, I), Z2) :-
    ExpR is exp(R),
    RotR is cos(I),
    RotI is sin(I),
    complex_mul(c(ExpR, 0), c(RotR, RotI), Z2).

:- inline(complex_abs/2).
complex_abs(c(R, I), Abs) :-
    Abs is sqrt(R*R + I*I).

complex_sqrt(c(R, 0), Z2) :-
    R >= 0, !, SqrtR is sqrt(R), Z2 = c(SqrtR, 0).

complex_sqrt(c(R, 0), Z2) :-
    R < 0, !, SqrtR is sqrt(-R), Z2 = c(0, SqrtR).

complex_sqrt(c(R, I), Z2) :-
    complex_abs(c(R, I), Abs),
    SqrtAbs is sqrt(Abs),
    atan2_raw(I, R, Td), T is Td/2,
    complex_exp(c(0, T), ExpT),
    complex_mul(c(SqrtAbs, 0), ExpT, Z2).

:- inline(complex_sinh/2).
complex_sinh(Z1, Z2) :-
    OneHalf = c(0.5, 0),
    complex_neg(Z1, NegZ1),
    complex_exp(NegZ1, ExpNegZ1),
    complex_neg(ExpNegZ1, NegExpNegZ1),
    complex_exp(Z1, ExpZ1),
    complex_add(ExpZ1, NegExpNegZ1, Twice),
    complex_mul(Twice, OneHalf, Z2).

:- inline(complex_sin/2).
complex_sin(Z1, Z2) :-
    I = c(0, 1),
    complex_mul(Z1, I, Sl1),
    complex_sinh(Sl1, Sl2),
    complex_mul(Sl2, I, Sl3),
    complex_neg(Sl3, Z2).

:- inline(complex_cosh/2).
complex_cosh(Z1, Z2) :-
    OneHalf = c(0.5, 0),
    complex_neg(Z1, NegZ1),
    complex_exp(NegZ1, ExpNegZ1),
    complex_exp(Z1, ExpZ1),
    complex_add(ExpZ1, ExpNegZ1, Twice),
    complex_mul(Twice, OneHalf, Z2).

:- inline(complex_cos/2).
complex_cos(Z1, Z2) :-
    I = c(0, 1),
    complex_mul(Z1, I, Sl),
    complex_cosh(Sl, Z2).

:- inline(complex_log/2).
complex_log(Z1, Z2) :-
    complex_log(Z1, Z2, 0).

:- inline(complex_log/3).
complex_log(Z1, Z2, K) :-
    c(R, I) = Z1,
    complex_abs(Z1, Rr),
    atan2_raw(I, R, T1),
    T is T1 + 2 * pi * K,
    LogR is ln(Rr),
    c(LogR, T) = Z2.

:- inline(complex_atanh/2).
complex_atanh(Z1, Z2) :-
    complex_atanh(Z1, Z2, 0).

:- inline(complex_atanh/3).
complex_atanh(Z1, Z2, K) :-
    complex_add(c(1, 0), Z1, Sl1),
    complex_add(c(-1, 0), Z1, Sl2),
    complex_neg(Sl2, Sl3),
    complex_recip(Sl3, Sl4),
    complex_mul(Sl1, Sl4, Res1),
    complex_log(Res1, Res2, K),
    complex_mul(c(0.5, 0), Res2, Z2).

:- inline(complex_atan/2).
complex_atan(Z1, Z2) :-
    complex_atan(Z1, Z2, 0).

:- inline(complex_atan/3).
complex_atan(Z1, Z2, K) :-
    complex_mul(Z1, c(0, 1), Sl),
    complex_atanh(Sl, Sl2, K),
    complex_mul(c(0, -1), Sl2, Z2).

:- inline(complex_tan/2).
complex_tan(Z1, Z2) :-
    complex_sin(Z1, SinZ1),
    complex_cos(Z1, CosZ1),
    complex_recip(CosZ1, CosZ1Recip),
    complex_mul(CosZ1Recip, SinZ1, Z2).

:- inline(complex_tanh/2).
complex_tanh(Z1, Z2) :-
    complex_sinh(Z1, SinhZ1),
    complex_cosh(Z1, CoshZ1),
    complex_recip(CoshZ1, CoshZ1Recip),
    complex_mul(CoshZ1Recip, SinhZ1, Z2).

:- inline(complex_asinh/2).
complex_asinh(Z1, Z2) :-
    complex_mul(Z1, Z1, Z1Square),
    complex_add(Z1Square, c(1, 0), Loc1),
    complex_sqrt(Loc1, Loc2),
    complex_add(Loc2, Z1, Loc3),
    complex_log(Loc3, Z2).

:- inline(complex_acosh/2).
complex_acosh(Z1, Z2) :-
    complex_mul(Z1, Z1, Z1Square),
    complex_add(Z1Square, c(-1, 0), Loc1),
    complex_sqrt(Loc1, Loc2),
    complex_add(Loc2, Z1, Loc3),
    complex_log(Loc3, Z2).

:- inline(complex_asin/2).
complex_asin(Z1, Z2) :-
    complex_mul(c(0, 1), Z1, Loc),
    complex_asinh(Loc, Loc1),
    complex_mul(c(0, -1), Loc1, Z2).

:- inline(complex_acos/2).
complex_acos(Z1, Z2) :-
    HalfPi is pi / 2,
    complex_asin(Z1, Loc),
    complex_neg(Loc, Loc1),
    complex_add(c(HalfPi, 0), Loc1, Z2).

atan2_raw(Y, X, Res) :-
        Res is atan(Y, X).

complex_sum([], c(0, 0)).
complex_sum([H], H).
complex_sum(Ss, c(Rt, It)) :-
        (foreach(R, Rs), foreach(I, Is), foreach(S, Ss)
        do S = c(R, I)),
        Rt is sum(Rs), It is sum(Is).

complex_prod([], c(1, 0)).
complex_prod(L, Res) :-
        complex_prod0(L, c(1, 0), c(R_expr, I_expr)),
        R is eval(R_expr),
        I is eval(I_expr),
        Res = c(R, I).

complex_prod0([], Res, Res).
complex_prod0([c(R0, I0) | T], c(Ra, Ia), Res) :-
        Rn = c(R0 * Ra - I0 * Ia, Ra * I0 + R0 * Ia),
        complex_prod0(T, Rn, Res).


complex_pow(X, Z, Y) :-
        Z = _ c _,
        complex_log(X, LogX),
        complex_mul(LogX, Z, ZLogX),
        complex_exp(ZLogX, Y).
complex_pow(_, 0, c(1, 0)).
complex_pow(X, 1, X).
complex_pow(X, N, Y) :-
        integer(N), N < 0,
        Xr c_is / X,
        NegN is -N,
        complex_pow(Xr, NegN, Y).
complex_pow(X, N, Y) :-
        integer(N), N >= 2,
        power_list(N, X, L),
        complex_prod(L, Y).
        

c_is(c(X1, X2), c(X3, X4)) :-
        X1 is X3,
        X2 is X4,
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
        Zneg c_is - Ze,
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

complex_sqr(X1 c X2, Y1 c Y2) :-
        Y1 is X1^2 - X2^2,
        Y2 is 2*X1*X2.

power_list(N, X, L) :-
        (fromto([N, X, L], In, Out, [0, _, []])
        do In = [N, Xi, L],
           Nd is N // 2,
           Nr is N rem 2,
           (Nr =:= 0 -> Ln = L; L = [Xi | Ln]),
           complex_sqr(Xi, Xip),
           Out = [Nd, Xip, Ln]).

binom(N, K, R) :-
        (for(I, 1, K, 1),
         fromto([1, 1], [In1, In2],
                [Out1, Out2], [Num, Den]), param(N)
        do  Out1 is In1 * (N - I + 1),
            Out2 is In2 * I),
        R is Num // Den.