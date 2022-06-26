import math

import sympy as sp
x = sp.symbols('x')

polinom = sp.sin(x**2+2)-sp.ln(2*x**5-3*x) #function
start_point = 2.5 * math.pi
end_point = 5.6 * math.pi
n = 16 #number of R(n,n) in Rumberg

def rumbergMethod(polinom,start_point,end_point,n):
    """
    Finds and returns the result of integral by Rumberg method
    :param polinom: function on which apply the integral
    :param start_point: start point
    :param end_point: end point
    :param n: number of R(n,n) in Rumberg
    :return: the result of integral by Rumberg method
    """
    f = sp.lambdify(x, polinom)
    i = 0
    RArray = []
    print("Building R Vector:")
    while i < n:
        multiplicity = (end_point-start_point)*(0.5**(i+1))
        sum = 0
        #calculate sum
        t = 0
        if i != 0:
            end_of_sum = (2 ** i) - 1
            sum_section_size = (end_point-start_point)/(end_of_sum+1)
            curr_point = start_point + sum_section_size
            while t < end_of_sum:
                sum += 2 * f(curr_point)
                curr_point += sum_section_size
                t += 1
        sum = sum + f(start_point) + f(end_point)
        #append result
        RArray.append(multiplicity * sum)
        print("R(" + str(i + 1) + ",1) = (" + str(end_point) + " - " +str(start_point) + ") * (0.5 ** " +str(i+1) +") * " + str(sum) + " = " + str(RArray[i]))
        i += 1
    print("-------")
    tempRArray = []
    j = 1
    while j < n:
        i = 1
        while i < n - j + 1:
            tempRArray.append(RArray[i] + (1/((4 ** j)-1) * (RArray[i] - RArray[i-1])))
            print("R(" + str(i + j) + "," + str(j + 1) + ") = R(" + str(i + j) + "," + str(j) + ") + (1/((4 ** "+str(j)+")-1) * (R(" + str(i + j) + "," + str(j) + ") - R(" + str(i + j - 1) + "," + str(j) + ")) = " + str(tempRArray[i-1]))
            i += 1
        j += 1
        RArray = tempRArray.copy()
        tempRArray.clear()
        print("-------")
    return RArray[0]

print("Integral by Rumberg method:")
print(rumbergMethod(polinom,start_point,end_point,n))