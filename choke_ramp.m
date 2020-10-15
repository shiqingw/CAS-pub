function ch_op = choke_ramp(current_pos,set_pos,start_t,duration,t,dt)

ramp = (set_pos - current_pos)/duration;

if t*dt<start_t
    ch_op = current_pos;
end

if (t*dt >= start_t) && (t*dt <= start_t+duration)
    if ramp > 0
        ch_op = max(current_pos, current_pos + ramp*(t*dt-start_t));
    else
        ch_op = max(set_pos, current_pos + ramp*(t*dt-start_t));
    end
end

if t*dt>start_t+duration
    ch_op = set_pos;
end

end
