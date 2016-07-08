function [u2] = mutation(u)
    x=randi(length(u));
    u2=u;
    if (u(x)==0)
        u2(x)=1;
    else
        u2(x)=0;
    end
end




