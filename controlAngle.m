function theta = controlAngle(yP,uconP,thetamax)

if yP<uconP
    theta=asin(yP/uconP); %control angle, pur
    if abs(theta)>=thetamax
        theta=sign(theta)*thetamax;
    end
else %if yP>uconP
    if yP>0
        theta=thetamax;
    else %if yP<0
        theta=-thetamax;
    end
end

end

