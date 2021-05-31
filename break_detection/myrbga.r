## Code provided by Andrea Toreti

myrbga<-function (stringMin = c(), stringMax = c(), suggestions = NULL, 
    popSize = 50, iters = 100, mutationChance = 0.6, elitism = NA, 
    monitorFunc = NULL, evalFunc = NULL, showSettings = FALSE, 
    verbose = FALSE,ydati,ypi) 
{


    vars<-length(stringMin)
    K<-stringMin[vars]





        elitism<-floor(popSize*0.1)


    if (vars > 0) {
        if (!is.null(suggestions)) {

            population = matrix(nrow = popSize, ncol = vars)
            suggestionCount = dim(suggestions)[1]
            for (i in 1:suggestionCount) {
                population[i, ] = suggestions[i, ]
            }

            for (var in 1:vars) {
                population[(suggestionCount + 1):popSize, var] = stringMin[var] + 
                  round(runif(popSize - suggestionCount),digit=0) * (stringMax[var] - 
                    stringMin[var])
            }
        }
        else {

            population = matrix(nrow = popSize, ncol = vars)

             mypop<-matrix(nrow=popSize,ncol=(K+1))
             mypop[,1]<-1
             mypop[,(K+1)]<-vars	    
             for(var in 1:popSize) {
                mypop[var,2:K]<-sort(sample(2:(vars-1),(K-1)))                                      
                 for(p in 1:K) {
                   population[var,mypop[var,p]:mypop[var,(p+1)]]<-p
                 } 
             }
        }

        bestEvals = rep(NA, iters)
        meanEvals = rep(NA, iters)
        evalVals = rep(NA, popSize)
        pevalVals<-rep(NA,popSize)
        for (iter in 1:iters) {


            if(iter==1){
             for (object in 1:popSize) {
                if (is.na(evalVals[object])) {

                  evalVals[object] = evalFunc(population[object,],a=ypi,b=ydati)

                }
              }
            } 
            bestEvals[iter] = min(evalVals)
            meanEvals[iter] = mean(evalVals)

            if (!is.null(monitorFunc)) {

                result = list(type = "floats chromosome", stringMin = stringMin, 
                  stringMax = stringMax, popSize = popSize, iter = iter, 
                  iters = iters, population = population, elitism = elitism, 
                  mutationChance = mutationChance, evaluations = evalVals, 
                  best = bestEvals, mean = meanEvals)
                class(result) = "rbga"
                monitorFunc(result)
            }
            if (iter < iters) {

                newPopulation = matrix(nrow = popSize, ncol = vars)
                newEvalVals = rep(NA, popSize)

                sortedEvaluations = sort(evalVals, index = TRUE)
                sortedPopulation = matrix(population[sortedEvaluations$ix, 
                  ], ncol = vars)


                intPopulation = matrix(nrow = popSize, ncol = vars)
                intEvalVals = rep(NA, popSize)

                for (tour in 1: popSize) {
                 champ<-sample(1:popSize,2)
                 if(evalVals[champ[1]]>= evalVals[champ[2]]){
                   intPopulation[tour,]<-population[champ[2],]
                 }
                 else {intPopulation[tour,]<-population[champ[1],]}
                } 


                if (vars > 1) {
                 child<-1
                  while (child <= popSize) {

                    parentIDs = sample(1:popSize, 2)
                    parents = intPopulation[parentIDs, ]

                   if(all(parents[1,]==parents[2,]))
                   { newPopulation[child, ] = parents[1, ]
                   child<-child+1
                   next
                   }
                    crossOverPoint = sample(1:(vars-1), 1)

                   crossProb<-runif(1)
                   if((crossProb<= 0.55)&((parents[2,(crossOverPoint+1)]-parents[1,crossOverPoint])<=1)&((parents[2,(crossOverPoint+1)]-parents[1,crossOverPoint])>=0)){
                      newPopulation[child, ] = c(parents[1, ][1:crossOverPoint], 
                        parents[2, ][(crossOverPoint + 1):vars])
                   child<-child+1
                   next                    
                   }
                   if((crossProb>0.55)|((parents[2,(crossOverPoint+1)]-parents[1,crossOverPoint])>1)|((parents[2,(crossOverPoint+1)]-parents[1,crossOverPoint])<0)){  
                     newPopulation[child, ]<-parents[1,]
                     if(child!=popSize) {
                     newPopulation[(child+1), ]<-parents[2,]
                     }
                    child<-child+2
                   }
                  }
                }



                if (mutationChance > 0) {

                  mutationCount = 0
                  for (object in 1:popSize) {
                    for (var in 2:(vars-1)) {
                     if(newPopulation[object,var]==(newPopulation[object,(var+1)]-1)){
                      if (runif(1) < mutationChance) {
                       newPopulation[object,var]<-newPopulation[object,var]+1
                      }
                     }
                    }
                  }


                }

             for (object in 1:popSize) {
                if (is.na(newEvalVals[object])) {
                  newEvalVals[object] = evalFunc(newPopulation[object,],a=ypi,b=ydati)
                }
             }                 

            bestnewEvals<-min(newEvalVals)
                 population<-newPopulation
                 evalVals<-newEvalVals

             if(bestnewEvals>bestEvals[iter]){ 

                newsortedEvaluations = sort(newEvalVals, index = TRUE)
                newsortedPopulation = matrix(newPopulation[newsortedEvaluations$ix, 
                  ], ncol = vars)
                 if (elitism > 0) {

                   newsortedPopulation[(popSize-elitism+1):popSize, ] = sortedPopulation[1:elitism, 
                     ]
                  newsortedEvaluations$x[(popSize-elitism+1):popSize] = sortedEvaluations$x[1:elitism]
                 }              

                 population = newsortedPopulation
                 evalVals = newsortedEvaluations$x

             }

            }
        }
    }
    result = list(type = "floats chromosome", stringMin = stringMin, 
        stringMax = stringMax, popSize = popSize, iters = iters, 
        suggestions = suggestions, population = population, elitism = elitism, 
        mutationChance = mutationChance, evaluations = evalVals, 
        best = bestEvals, mean = meanEvals)
    class(result) = "rbga"

    return(result)
}
