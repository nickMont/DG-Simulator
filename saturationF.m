function v = saturationF(v_in,vmin,vmax)

if v_in>=vmax
    v=vmax;
elseif v_in<=vmin
    v=vmin;
else
    v=v_in;
end

if vmin>=vmax
    fprintf('Error in saturation \n') %called saturation with max as min
end
        
end

