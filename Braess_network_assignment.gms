$title traffic assginment for Braess network

set i node /1*4/
    a link /1*5/
    p path /1*3/;

alias(i,j);

parameter
      q(i,j) the demand for od pair ij
             /1.2   4000/
      path_od(i,j,p) the pth path for od pair ij
             /1.2.1 1
              1.2.2 1
              1.2.3 1/
      fie(i,j,p,a) the pth path on od pair ij passes link a
             /1.2.1.1 1
              1.2.1.4 1
              1.2.2.1 1
              1.2.2.5 1
              1.2.2.2 1
              1.2.3.3 1
              1.2.3.2 1/
       cap(a) the capacity on link a
             /1 3000
              2 2000
              3 5000
              4 2000
              5 1000/
       FFTT(a) the free-flow travel time of lnk a
             /1 30
              2 45
              3 45
              4 20
              5 10/ ;

scalar beta the paremeter in BPR function  /4/
       alpha the parameter in BPR function /0.15/;

variable c(i,j,p) the travle time on path a for od pair ij
         TT(a) travel time on link a
         gap  the gap between all paths and mimimum travel time path for all od pair;

positive variable  link_flow(a)   the flow on link
                   r(i,j,p)       the flow on path p for od pair ij
                   pie(i,j)       the minima travel time for od pair ij;

equations
       Travl_Time(a) the travel time on link a
       od_flow_constraint(i,j)  the sum of all path flow on od pair ij equals to the od demand
       link_flow_constraint(a) the sum of all path flow passing link a equals the flow on link a
       path_traveltime_constraint(i,j,p) the travel time on path p for od pair ij equlas the sum of travel time on links passed by this path
       capacity_constraint(a) the flow on link a can't larger than the capacity
       pie_constraint(i,j,p)  the minima travel time for od pair ij
       obj  the objective function       ;
Travl_Time(a)..TT(a)=e=FFTT(a)*(1+alpha*(link_flow(a)/cap(a))**beta);
od_flow_constraint(i,j)$(q(i,j)>0)..sum(p,r(i,j,p))=e=q(i,j);
link_flow_constraint(a)..link_flow(a)=e=sum((i,j,p),r(i,j,p)*fie(i,j,p,a));
path_traveltime_constraint(i,j,p)$(path_od(i,j,p)=1)..c(i,j,p)=e=sum(a,TT(a)*fie(i,j,p,a));
capacity_constraint(a)..link_flow(a)=l=cap(a);
pie_constraint(i,j,p)$(q(i,j)>0)..pie(i,j)=l=c(i,j,p);
obj..gap=e=sum((i,j,p),r(i,j,p)*(c(i,j,p)-pie(i,j)));

model TABN /all/;
solve TABN using nlp minimizing gap;
display pie.l,link_flow.l,TT.l,r.l,c.l,gap.l;







