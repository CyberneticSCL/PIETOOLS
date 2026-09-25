function st = bl_settings(name)
st = lpisettings(name);
st.sos_opts.solver   = 'mosek';
st.sos_opts.simplify = false;
end


