# lb_compare.awk, 10/28/2026: stock vs sopvar, one line per (example,
# settings file) present in both result files.
# Usage: awk -F, -f lb_compare.awk catalog.csv results.csv results_sop.csv
# Output: comparison CSV on stdout; summary lines prefixed '#' on stderr.
# Certified: certified, converged, hi_pass, max_solves (a bisection that
# certified its lower bracket). Better result: a larger certified parameter
# (bisection), a smaller bound (gain or norm).
function cert(s) { return s=="certified" || s=="converged" || s=="hi_pass" || s=="max_solves" }
FILENAME ~ /catalog/ { if (FNR>1) { typ[$1"|"$2]=$4; bis[$1"|"$2]=$8; dim[$1"|"$2]=$11 } ; next }
FILENAME ~ /results_sop/ { if (FNR>1) { k=$1"|"$2"|"$3; ss[k]=$4; sr[k]=$5; st[k]=$8; sd[k]=$9; sn[k]=$10; doc[k]=$6 } ; next }
FILENAME ~ /results/ { if (FNR>1) { k=$1"|"$2"|"$3; os[k]=$4; or_[k]=$5; ot[k]=$8; od[k]=$9 } ; next }
END {
  print "library,id,settings,type,bisect,dim,documented,stock_status,stock_result,stock_seconds,stock_decvars,sop_status,sop_result,sop_seconds,sop_decvars,verdict"
  for (k in ss) {
    if (!(k in os)) continue
    split(k,a,"|"); c=a[1]"|"a[2]; t=typ[c]; b=bis[c]
    v=""
    if (ss[k]=="error" && os[k]=="error") v="both error"
    else if (ss[k]=="error") v="sop error"
    else if (os[k]=="error") v="stock error"
    else if (ss[k]=="skipped" || os[k]=="skipped") v="skipped"
    else if (cert(ss[k]) && !cert(os[k])) v="sop only"
    else if (!cert(ss[k]) && cert(os[k])) v="stock only"
    else if (!cert(ss[k]) && !cert(os[k])) v="neither"
    else {
      x=sr[k]+0; y=or_[k]+0; rel=(x-y)/((y<0?-y:y)+1e-12)
      if (t ~ /^stability/ && b!="yes") v="same"
      else if (rel<1e-6 && rel>-1e-6) v="same"
      else if (b=="yes") v=(x>y ? "sop tighter" : "sop looser")
      else v=(x<y ? "sop tighter" : "sop looser")
    }
    n[a[3]" "v]++
    if (ss[k]!="skipped" && os[k]!="skipped") { tsop[a[3]]+=st[k]; tst[a[3]]+=ot[k]; if ((sd[k]+0)>0 && (od[k]+0)>0) { dr[a[3]]+=log(sd[k]/od[k]); nd[a[3]]++ } }
    print a[1]","a[2]","a[3]","t","b","dim[c]","doc[k]","os[k]","or_[k]","ot[k]","od[k]","ss[k]","sr[k]","st[k]","sd[k]","v
  }
  for (q in n) print "# " q ": " n[q] > "/dev/stderr"
  for (q in tsop) printf "# %s time stock %.0f s, sop %.0f s; decision variables sop/stock geometric mean %.2f over %d rows\n", q, tst[q], tsop[q], exp(dr[q]/nd[q]), nd[q] > "/dev/stderr"
}
