##### Design different exp. design for using coef
## compare three groups
cond<- factor(rep(c("A","B", "C"), each=2))
model.matrix(~ cond)

## to compare B vs C, so B needs to be ref. level
cond<- relevel(cond, "B")
model.matrix(~ cond)

## Two factors and interaction
group<- factor(rep(1:3, each=4))
cond2<- factor(rep(rep(c("A", "B"), each=2),3))
model.matrix(~ group + cond2 + group:cond2)

## To compare between group3 vs group2,group2 needs to be made ref level; no change in condition
group<- relevel(group, "2")
model.matrix(~ group + cond2 + group:cond2)

## Two groups, two individuals per group, compare within-individual condition effects
group<- factor(rep(1:2,each=4))
ind<- factor(rep(rep(1:2,each=2),2))
cond<- factor(rep(c("A","B"),4))
model.matrix(~ group + group:ind + group:cond)

## To compare cond effect across group; cond would be a main effect
model.matrix(~ group + cond + group:ind + group:cond)
