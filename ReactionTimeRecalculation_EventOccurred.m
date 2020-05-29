function reactiontime = ReactionTimeRecalculation_EventOccurred( x, y, z, subvolume_nH, k, t_H)
    
    if subvolume_nH(x, y, z) == 0
        reactiontime = inf;
    else
        reactiontime = log(1/rand)/(6*subvolume_nH(x, y, z)*k) + t_H;
    end
    
end