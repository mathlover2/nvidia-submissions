:- module(complex_polar).
:- lib(ic).
:- local struct(c(re, im)).
:- local struct(cp(mag, arg)).

:- export re/2, im/2, complex_add/3, complex_mul/3.
:- export complex_neg/2, polar_form/2, complex_conjugate/2.
:- export complex_exp/2, complex_abs/2, complex_recip/2.
:- export complex_sqrt/2, complex_sinh/2, complex_sin/2.
:- export complex_cosh/2, complex_cos/2, complex_log/2, complex_log/3.
:- export complex_atanh/2, complex_atan/2, complex_tan/2.
:- export complex_tanh/2, complex_asinh/2, complex_acosh/2.
:- export complex_asin/2, complex_acos/2, atan2_raw/3.
:- export complex_eq/2, complex_sum/2, complex_prod/2.

re(c(R, _), R) :- !.
re(cp(Mag, Arg), R) :-
        R $= Mag * cos(Arg), !.

im(c(_, I), I) :- !.
im(cp(Mag, Arg), I) :-
        I $= Mag * sin(Arg), !.

polar_form(c(R, I), cp(Mag, Arg)) :-
        complex_abs(c(R, I), Mag),
        complex_exp(c(0, Arg), Swing),
        complex_mul(c(Mag,0), Swing, Test),
        complex_eq(c(R, I), Test).
        
complex_eq(c(R1, R2), c(I1, I2)) :-
        R1 $= I1,
        R2 $= I2, !.

complex_eq(cp(Mag1, _), cp(Mag2,_)) :-
        Mag1 $= 0,
        Mag2 $= 0, !.

complex_eq(cp(Mag1, Arg1), cp(Mag2, Arg2))  :-
        Mag1 $= Mag2,
        integers([K]),
        Arg1 - Arg2 $= 2*pi*K.

complex_eq(cp(Mag1, Arg1), cp(Mag2, Arg2))  :-
        Mag1 $= -Mag2,
        integers([K]),
        Arg1 - Arg2 $= (2*pi*K + pi).

r_to_c(R, c(R, 0)).
r_to_cp(R, cp(R, 0)) :- R $> 0.
r_to_cp(R, cp(Rn, Pi)) :- 
        R $< 0,
        Rn $= -R,
        Pi $= pi.

complex_add(c(R1, I1), c(R2, I2), c(R3, I3)) :-
    R3 $= R1 + R2,
    I3 $= I1 + I2.

complex_add(cp(M1, A1), cp(M2, A2), cp(M3, A3)) :-
        polar_form(X1, cp(M1, A1)),
        polar_form(X2, cp(M2, A2)),
        complex_add(X1, X2, X3),
        polar_form(X3, cp(M3, A3)).

complex_sum([], c(0, 0)).
complex_sum([H], H).
complex_sum(Ss, c(Rt, It)) :-
        (foreach(R, Rs), foreach(I, Is), foreach(S, Ss)
        do S = c(R, I)),
        Rt $= sum(Rs), It $= sum(Is).

complex_mul(c(R1, I1), c(R2, I2), c(R3, I3)) :-
    R3 $= R1 * R2 - I1 * I2,
    I3 $= R2 * I1 + R1 * I2.

complex_prod([], c(1, 0)).
complex_prod(L, Res) :-
        complex_prod0(L, c(1, 0), c(R_expr, I_expr)),
        R $= eval(R_expr),
        I $= eval(I_expr),
        Res = c(R, I).

complex_prod0([], Res, Res).
complex_prod0([c(R0, I0) | T], c(Ra, Ia), Res) :-
        Rn = c(R0 * Ra - I0 * Ia, Ra * I0 + R0 * Ia),
        complex_prod0(T, Rn, Res).

complex_neg(c(R1, I1), c(R2, I2)) :-
    R2 $= -R1,
    I2 $= -I1, !.

complex_conjugate(c(R1, I1), c(R1, I2)) :-
    I2 $= -I1.

complex_recip(c(R1, I1), c(R2, I2)) :-
    Norm $= R1*R1 + I1*I1,
    R2 $= R1 / Norm,
    I2 $= -I1 / Norm, !.

complex_exp(c(R, I), Z2) :-
    ExpR $= exp(R),
    RotR $= cos(I),
    RotI $= sin(I),
    complex_mul(c(ExpR, 0), c(RotR, RotI), Z2).

complex_abs(c(R, I), Abs) :-
    Abs $= sqrt(R*R + I*I), !.

complex_sqrt(c(R, 0), Z2) :-
    R $>= 0, SqrtR $= sqrt(R), Z2 = c(SqrtR, 0), !.

complex_sqrt(c(R, 0), Z2) :-
    R $< 0, SqrtR $= sqrt(-R), Z2 = c(0, SqrtR), !.

complex_sqrt(c(R, I), Z2) :-
    complex_abs(c(R, I), Abs),
    SqrtAbs $= sqrt(Abs),
    atan2_raw(I, R, Td), T $= Td/2,
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

complex_log(Z1, Z2, K) :-
    c(R, I) = Z1,
    complex_abs(Z1, Rr),
    atan2_raw(I, R, T1),
    T $= T1 + 2 * pi * K,
    LogR $= ln(Rr),
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
    complex_mul(Z1, Z1, Z1Square),
    complex_add(Z1Square, c(1, 0), Loc1),
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
    HalfPi $= pi / 2,
    complex_asin(Z1, Loc),
    complex_neg(Loc, Loc1),
    complex_add(c(HalfPi, 0), Loc1, Z2).

% atan2_raw(Y, X, Res) :-
%         X $> 0, Res $= atan(Y/X).
% atan2_raw(Y, X, Res) :-
%         X $< 0, Y $>= 0, Res $= atan(Y/X) + pi.
% atan2_raw(Y, X, Res) :-
%         X $< 0, Y $< 0, Res $= atan(Y/X) - pi.
% atan2_raw(Y, X, Res) :-
%         X $= 0, Y $> 0, Res $= pi/2.
% atan2_raw(Y, X, Res) :-
%         X $= 0, Y $< 0, Res $= -pi/2.

atan2_raw(Y, X, Res) :-
        (X $> 0 and Res $= atan(Y/X))
        or (X $< 0 and Y $>= 0 and Res $= atan(Y/X) + pi)
        or (X $< 0 and Y $< 0 and Res $= atan(Y/X) - pi)
        or (X $= 0 and Y $> 0 and Res $= pi/2)
        or (X $= 0 and Y $< 0 and Res $= -pi/2).

atan2_raw(Y, X, _) :-
        X $= 0 and Y $= 0, fail.
