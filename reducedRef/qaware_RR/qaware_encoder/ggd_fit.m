function best_paras = ggd_fit(hist, width)

best_paras = abs(fminsearch('ggd_error', [2 0.6], [], hist, width));

return
