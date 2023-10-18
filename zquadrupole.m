function Gz=zquadrupole(z,G,b,L,s)
function gz=quadfunshift(z,G,b,L,s) %smoothfield %s as a shift parameter to recenter field
    gz=G/2*(tanh(1/2*b*(1/2*L-(z-s))))+G/2*(tanh(1/2*b*(1/2*L+(z-s))));
end
function gz=hardquadfunshift(z,G,L,s) %hardedge
    gz=G*(heaviside(z-s+L/2)-heaviside(z-s-L/2));
end
if b~=0
    Gz=quadfunshift(z,G,b,L,s);
else
    Gz=hardquadfunshift(z,G,L,s);
end

end