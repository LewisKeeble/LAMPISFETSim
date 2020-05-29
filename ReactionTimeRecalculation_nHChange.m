function reactiontime_out = ReactionTimeRecalculation_nHChange( x, y, z, subvolume_nH, reactiontime_in, k_diffusion, t_H)
    
    if subvolume_nH(x, y, z) == 0
        reactiontime_out = inf;
    elseif subvolume_nH(x, y, z) == 1
        reactiontime_out = log(1/rand)/(6*k_diffusion) + t_H;
    else
        reactiontime_out = ((subvolume_nH(x, y, z)-1)/subvolume_nH(x, y, z))*(reactiontime_in - t_H) + t_H;
    end
                    
end