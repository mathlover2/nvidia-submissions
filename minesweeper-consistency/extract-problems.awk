function make_substitutions(string,    f) {
    f = gensub(/_/, "-2", "g", string);
    f = gensub(/]],/, "|];", "g", f);
    f = gensub(/],/, ",", "g", f);
    return f;
}

BEGIN {
    OFS="";
    printing_flag = 0;
    problem_number = 1;
    currfile = "";
}

/^problem\(/ {
    printing_flag = 1;
    currfile = "minesweeper_consistency_" problem_number ".dzn";
    f = gensub(/\[\[/, "[|", "g", $2);
    print "initial=", make_substitutions(f) >> currfile;
    next;
}



/\[[[:digit:]]+, [[:digit:]]+\]\)\./ {
    if (printing_flag == 1) {
	f = gensub(/\[/, "|", "g", $1);
	print make_substitutions(f) >> currfile;
	match($0, /\[([[:digit:]]+), ([[:digit:]]+)\]/, numpairstring);
	print "h=", numpairstring[1], ";" >> currfile;
	print "w=", numpairstring[2], ";" >> currfile;
	problem_number++;
	printing_flag = 0;
	next;
    }
}

( printing_flag == 1 ) {
    f = gensub(/\[/, "|", "g", $1);
    print make_substitutions(f) >> currfile; 
}



