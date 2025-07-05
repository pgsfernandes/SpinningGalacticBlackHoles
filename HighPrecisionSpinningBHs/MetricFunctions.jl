function gtt(x, y, f::Field, g::Field, h::Field, W::Field,rh::Double64=1.0;dx::Integer=0,dr::Bool=false)
    if dx==0
        return -(((1 + x)^2*f(x,y))/(-3 + x)^2) + ((-3 + x)^4*(-1 + x)^2*g(x,y)*sin(y)^2*W(x,y)^2)/(64*rh^2*f(x,y))
    elseif dx==1
        if !dr
            return (-128*(-3 + x)*(1 + x)*f(x,y) + 128*(1 + x)^2*f(x,y) + (2*(-3 + x)^7*(-1 + x)*g(x,y)*sin(y)^2*W(x,y)^2)/(rh^2*f(x,y)) + (4*(-3 + x)^6*(-1 + x)^2*g(x,y)*sin(y)^2*W(x,y)^2)/(rh^2*f(x,y)) - 64*(-3 + x)*(1 + x)^2*f(x,y,dx=1) - ((-3 + x)^7*(-1 + x)^2*g(x,y)*sin(y)^2*W(x,y)^2*f(x,y,dx=1))/(rh^2*f(x,y)^2) + ((-3 + x)^7*(-1 + x)^2*sin(y)^2*W(x,y)^2*g(x,y,dx=1))/(rh^2*f(x,y)) + (2*(-3 + x)^7*(-1 + x)^2*g(x,y)*sin(y)^2*W(x,y)*W(x,y,dx=1))/(rh^2*f(x,y)))/(64*(-3 + x)^3)
        else
            return (1-x)^2/(2*rh) * gtt(x,y,f,g,h,W,rh,dx=1,dr=false)
        end
    else
        if !dr
            return (-384*(1 + x)^2*f(x,y) + (2*(-3 + x)^6*(39 - 50*x + 15*x^2)*g(x,y)*sin(y)^2*W(x,y)^2)/(rh^2*f(x,y)) + 256*(-3 + x)*(1 + x)*(2*f(x,y) + (1 + x)*f(x,y,dx=1)) + (4*(-3 + x)^7*(-1 + x)*(-5 + 3*x)*sin(y)^2*W(x,y)*(f(x,y)*W(x,y)*g(x,y,dx=1) + g(x,y)*(-(W(x,y)*f(x,y,dx=1)) + 2*f(x,y)*W(x,y,dx=1))))/(rh^2*f(x,y)^2) - 64*(-3 + x)^2*(2*f(x,y) + (1 + x)*(4*f(x,y,dx=1) + (1 + x)*f(x,y,dx=2))) + ((-3 + x)^8*(-1 + x)^2*sin(y)^2*(f(x,y)*W(x,y)*(4*f(x,y)*g(x,y,dx=1)*W(x,y,dx=1) + W(x,y)*(-2*f(x,y,dx=1)*g(x,y,dx=1) + f(x,y)*g(x,y,dx=2))) + g(x,y)*(2*f(x,y)^2*W(x,y,dx=1)^2 + W(x,y)^2*(2*f(x,y,dx=1)^2 - f(x,y)*f(x,y,dx=2)) + 2*f(x,y)*W(x,y)*(-2*f(x,y,dx=1)*W(x,y,dx=1) + f(x,y)*W(x,y,dx=2)))))/(rh^2*f(x,y)^3))/(64*(-3 + x)^4)
        else
            return (-1+x)^3/(4*rh^2) * ( 2*gtt(x,y,f,g,h,W,rh,dx=1,dr=false) + (-1+x)*gtt(x,y,f,g,h,W,rh,dx=2,dr=false) )
        end
    end
end

function gtphi(x, y, f::Field, g::Field, h::Field, W::Field,rh::Double64=1.0;dx::Integer=0,dr::Bool=false)
    if dx==0
        return -1/16*((-3 + x)^4*g(x,y)*sin(y)^2*W(x,y))/f(x,y)
    elseif dx==1
        if !dr
            return -1/16*((-3 + x)^3*sin(y)^2*(-((-3 + x)*g(x,y)*W(x,y)*f(x,y,dx=1)) + f(x,y)*((-3 + x)*W(x,y)*g(x,y,dx=1) + g(x,y)*(4*W(x,y) + (-3 + x)*W(x,y,dx=1)))))/f(x,y)^2
        else
            return (1-x)^2/(2*rh) * gtphi(x,y,f,g,h,W,rh,dx=1,dr=false)
        end
    else
        if !dr
            return ((-3 + x)^2*sin(y)^2*(2*(-3 + x)*f(x,y)*(-4*f(x,y) + (-3 + x)*f(x,y,dx=1))*(W(x,y)*g(x,y,dx=1) + g(x,y)*W(x,y,dx=1)) - g(x,y)*W(x,y)*(12*f(x,y)^2 + 2*(-3 + x)^2*f(x,y,dx=1)^2 - (-3 + x)*f(x,y)*(8*f(x,y,dx=1) + (-3 + x)*f(x,y,dx=2))) - (-3 + x)^2*f(x,y)^2*(2*g(x,y,dx=1)*W(x,y,dx=1) + W(x,y)*g(x,y,dx=2) + g(x,y)*W(x,y,dx=2))))/(16*f(x,y)^3)
        else
            return (-1+x)^3/(4*rh^2) * ( 2*gtphi(x,y,f,g,h,W,rh,dx=1,dr=false) + (-1+x)*gtphi(x,y,f,g,h,W,rh,dx=2,dr=false) )
        end
    end
end

function gphiphi(x, y, f::Field, g::Field, h::Field, W::Field,rh::Double64=1.0;dx::Integer=0,dr::Bool=false)
    if dx==0
        return (rh^2*(-3 + x)^4*g(x,y)*sin(y)^2)/(4*(-1 + x)^2*f(x,y))
    elseif dx==1
        if !dr
            return (rh^2*(-3 + x)^3*sin(y)^2*(-((3 - 4*x + x^2)*g(x,y)*f(x,y,dx=1)) + f(x,y)*(2*(1 + x)*g(x,y) + (3 - 4*x + x^2)*g(x,y,dx=1))))/(4*(-1 + x)^3*f(x,y)^2)
        else
            return (1-x)^2/(2*rh) * gphiphi(x,y,f,g,h,W,rh,dx=1,dr=false)
        end
    else
        if !dr
            return (rh^2*(-3 + x)^2*sin(y)^2*(2*(9 + 2*x + x^2)*f(x,y)^2*g(x,y) + 4*(-3 + x)*(-1 + x)*(1 + x)*f(x,y)*(-(g(x,y)*f(x,y,dx=1)) + f(x,y)*g(x,y,dx=1)) + (-3 + x)^2*(-1 + x)^2*(g(x,y)*(2*f(x,y,dx=1)^2 - f(x,y)*f(x,y,dx=2)) + f(x,y)*(-2*f(x,y,dx=1)*g(x,y,dx=1) + f(x,y)*g(x,y,dx=2)))))/(4*(-1 + x)^4*f(x,y)^3)
        else
            return (-1+x)^3/(4*rh^2) * ( 2*gphiphi(x,y,f,g,h,W,rh,dx=1,dr=false) + (-1+x)*gphiphi(x,y,f,g,h,W,rh,dx=2,dr=false) )
        end
    end
end

function grr(x, y, f::Field, g::Field, h::Field, W::Field,rh::Double64=1.0;dx::Integer=0,dr::Bool=false)
    if dx==0
        return ((-3 + x)^4*g(x,y)*h(x,y))/(16*f(x,y))
    elseif dx==1
        if !dr
            return ((-3 + x)^3*(-((-3 + x)*g(x,y)*h(x,y)*f(x,y,dx=1)) + f(x,y)*((-3 + x)*h(x,y)*g(x,y,dx=1) + g(x,y)*(4*h(x,y) + (-3 + x)*h(x,y,dx=1)))))/(16*f(x,y)^2)
        else
            return (1-x)^2/(2*rh) * grr(x,y,f,g,h,W,rh,dx=1,dr=false)
        end
    else
        if !dr
            return (-2*(3 - x)^3*f(x,y)*(4*f(x,y) - (-3 + x)*f(x,y,dx=1))*(h(x,y)*g(x,y,dx=1) + g(x,y)*h(x,y,dx=1)) + (-3 + x)^2*g(x,y)*h(x,y)*(12*f(x,y)^2 + 2*(-3 + x)^2*f(x,y,dx=1)^2 - (-3 + x)*f(x,y)*(8*f(x,y,dx=1) + (-3 + x)*f(x,y,dx=2))) + (-3 + x)^4*f(x,y)^2*(2*g(x,y,dx=1)*h(x,y,dx=1) + h(x,y)*g(x,y,dx=2) + g(x,y)*h(x,y,dx=2)))/(16*f(x,y)^3)
        else
            return (-1+x)^3/(4*rh^2) * ( 2*grr(x,y,f,g,h,W,rh,dx=1,dr=false) + (-1+x)*grr(x,y,f,g,h,W,rh,dx=2,dr=false) )
        end
    end
end

function gthetatheta(x, y, f::Field, g::Field, h::Field, W::Field,rh::Double64=1.0;dx::Integer=0,dr::Bool=false)
    if dx==0
        return (rh^2*(-3 + x)^4*g(x,y)*h(x,y))/(4*(-1 + x)^2*f(x,y))
    elseif dx==1
        if !dr
            return (rh^2*(-3 + x)^3*(-((3 - 4*x + x^2)*g(x,y)*h(x,y)*f(x,y,dx=1)) + f(x,y)*((3 - 4*x + x^2)*h(x,y)*g(x,y,dx=1) + g(x,y)*(2*(1 + x)*h(x,y) + (3 - 4*x + x^2)*h(x,y,dx=1)))))/(4*(-1 + x)^3*f(x,y)^2)
        else
            return (1-x)^2/(2*rh) * gthetatheta(x,y,f,g,h,W,rh,dx=1,dr=false)
        end
    else
        if !dr
            return (rh^2*(-3 + x)^2*(2*(9 + 2*x + x^2)*f(x,y)^2*g(x,y)*h(x,y) + 4*(-3 + x)*(-1 + x)*(1 + x)*f(x,y)*(f(x,y)*h(x,y)*g(x,y,dx=1) + g(x,y)*(-(h(x,y)*f(x,y,dx=1)) + f(x,y)*h(x,y,dx=1))) - (-3 + x)^2*(-1 + x)^2*(2*f(x,y)*(g(x,y)*f(x,y,dx=1) - f(x,y)*g(x,y,dx=1))*h(x,y,dx=1) - h(x,y)*(g(x,y)*(2*f(x,y,dx=1)^2 - f(x,y)*f(x,y,dx=2)) + f(x,y)*(-2*f(x,y,dx=1)*g(x,y,dx=1) + f(x,y)*g(x,y,dx=2))) - f(x,y)^2*g(x,y)*h(x,y,dx=2))))/(4*(-1 + x)^4*f(x,y)^3)
        else
            return (-1+x)^3/(4*rh^2) * ( 2*gthetatheta(x,y,f,g,h,W,rh,dx=1,dr=false) + (-1+x)*gthetatheta(x,y,f,g,h,W,rh,dx=2,dr=false) )
        end
    end
end