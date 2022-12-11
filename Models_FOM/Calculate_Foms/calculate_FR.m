function FR=calculate_FR(elementdata,refdata)
FR=0.5*sum(abs(elementdata-refdata));
end