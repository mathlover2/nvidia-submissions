:- use_module(complex).
:- lib(ic).
:- local struct(range(lo, hi)).

exclude_range_search(Prec, Goal, VarsorCVars, Res) :-
        process_vars(VarsorCVars, Vars),
        bag_create(Data), bag_create(Sols),
        ((call(Goal),
          bag_retrieve(Data, ExclusionConstraintData),
          create_constraint(Vars, ExclusionConstraintData,
                            ExclusionConstraint),
          call(ExclusionConstraint),
          bag_enter(Sols, Vars),
          maplist(range_from_data(Prec), Vars, NewData),
          bag_enter(Data, NewData),
          fail); true),
        bag_abolish(Data),
        bag_dissolve(Sols, Res).

process_vars([], []).
process_vars([X1 c X2 | T1], [X1, X2 | T2]) :-
        process_vars(T1, T2), !.
process_vars([X | T1], [X | T2]) :-
        process_vars(T1, T2), !.

range_from_data(Prec, Var, range(Lo, Hi)) :-
        (ic:is_solver_var(Var) ->
            (ic:get_bounds(Var, L, H),
             Lo is float(L - Prec),
             Hi is float(H + Prec));
            (breal_bounds(Var, Min, Max),
             Lo is float(Min - Prec),
             Hi is float(Max + Prec))).


create_constraint(_, [], 1 $= 1).
create_constraint(Vars, [D0 | Dt], Con0 and Rest) :-
        constraint_from_ranges(D0, Vars, Con0),
        create_constraint(Vars, Dt, Rest).

constraint_from_ranges([], [], 1 $= 0).
constraint_from_ranges([range(Lo, Hi) | T], [V | Vt], (V $> Hi or V $< Lo) or Rest) :-
        constraint_from_ranges(T, Vt, Rest).
        
% This function takes a list of numerical vectors and deletes those
% "duplicates" in which each entry is "close" to a previously
% encountered entry, according to a function `trim_sols` which takes
% `Prec` as a parameter.
%
% This is of course not well-defined in general, but when dealing with
% isolated clusters of intervals where members of each cluster are
% overlapping, it works quite well.
trim_sols(Bag1, Prec, Bag) :-
        trim_sols(Prec, Bag1, [], Bag), !.

trim_sols(_, [], New, New).
trim_sols(Prec, [H | T], [], New) :-
        trim_sols(Prec, T, [H], New).
trim_sols(Prec, [H | T], Prev, New) :-
        (near_member(Prec, H, Prev)
        -> trim_sols(Prec, T, Prev, New) ;
           trim_sols(Prec, T, [H | Prev], New)).

near_member(Prec, H, [L]) :- near(Prec, H, L).
near_member(Prec, H, [L0 | T]) :-
        near(Prec, H, L0);
        near_member(Prec, H, T).

near(_, [], []).
near(Prec, [TestVar | TestVarT], [OldVar | OldVarT]) :-
        ic:get_bounds(OldVar, OldVarL, OldVarH),
        ic:get_bounds(TestVar, TestVarL, TestVarH),
        (abs(TestVarL - OldVarL) $=< Prec) 
        or (abs(TestVarH - OldVarH) $=< Prec),
        near(Prec, TestVarT, OldVarT).