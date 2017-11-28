#install packages if required
install.packages("ReinforcementLearning")
install.packages("fTrading")
install.packages("TTR")
install.packages("devtools")
install.packages("reticulate") 
devtools::install_github("rstudio/keras")
library(keras)
library(tensorflow)
#install_tensorflow(conda="tf=keras")
library(ReinforcementLearning)
library(fTrading)
library(TTR)
library(reticulate)

setwd("Z:/EMJMCMC2016/examples/reinforce bot")

paramsim = list(seed = 1, n = 10, n.curs = 500,buypercent=0.9)
parambot = list(limit_trade_pairs  =50, min_buy_size = 0.001, ban_days_new_pair = 10,days_to_keep = 10, ban_days = 7, min_profit = 0.02, sell_loss = -0.35, min_currency_price = 0.000001,  trade_step = 999, balance = 10000, x.hist = lapply(FUN = function(x) {y = cumsum(sample(c(-10, 10),size = x, TRUE));rands = abs(rnorm(n = x,mean = 0,sd = 1));return(list(High= y + abs(min(y))+rands,Low =  y + abs(min(y))-rands, vol.high = runif(1,0,1000),vol.low = runif(1,0,1000) ))},X = rep(list(x = paramsim$n),paramsim$n.curs)))
#simulate simplified bot trading (to be replaced with a real bot on real markets at the end)
simulateAction=function(paramsim,parambot)
{
  balance = parambot$balance
  #simulate paramsim$n.curs markets for paramsim$n time slots
  set.seed(paramsim$seed)
  x = lapply(FUN = function(x) {y = cumsum(c(parambot$x.hist$high,sample(c(-10, 10),size = x, TRUE)));rands = abs(rnorm(n = x,mean = 0,sd = 1));return(list(High= y + abs(min(y))+rands,Low =  y + abs(min(y))-rands, vol.high = runif(1,0,1000),vol.low = runif(1,0,1000) ))},X = rep(list(x = paramsim$n),paramsim$n.curs))
  #calculate volatilities for all of the pairs
  vols = unlist(lapply(FUN = function(x) {TTR::chaikinVolatility(HL = cbind(x$High,x$Low))[paramsim$n]},X = parambot$x.hist))
  curs.buy.ids = which(vols>=sort(vols,decreasing = T)[parambot$limit_trade_pairs])
  #simulate buying the coins with highest volatilities with the lowest prices
  buys<-array(0,c(parambot$limit_trade_pairs,2))
  for(i in 1:parambot$limit_trade_pairs)
  {
    if(parambot$x.hist[[i]]$Low[paramsim$n]>=parambot$min_currency_price && parambot$x.hist[[i]]$vol.low>=parambot$min_buy_size)
    { 
      buys[i,1] = parambot$x.hist[[i]]$vol.low
      buys[i,2] = parambot$x.hist[[i]]$Low[paramsim$n]
      balance = balance -  buys[i,1]* buys[i,2]
    }
  }
  #simulate selling the coins
  for(i in 1:parambot$limit_trade_pairs)
  {
    for(j in seq(1,paramsim$n,by = parambot$trade_step))
    {
      
      if((x[[i]]$High[j]- buys[i,2])/(buys[i,2])>=parambot$min_profit)
      {
        vol = min(x[[i]]$vol.high*paramsim$buypercent,buys[i,1])
        balance = balance + vol*x[[i]]$High[j]
        buys[i,1] =  buys[i,1] - vol
        if(buys[i,1]==0)
          break
      }else if((x[[i]]$Low[j]- buys[i,2])/(buys[i,2])<=parambot$sell_loss)
      {
        balance = balance + buys[i,1]*x[[i]]$Low[j]
        buys[i,1] = 0
        buys[i,2] = 0
        break
      }
    }
  }
  balance = balance + sum(buys[,1]* buys[,2])
  return(list(balance=balance,profit = (balance-parambot$balance)/abs(parambot$balance),x = x))
}

#simulate some initial data
N=15001
balance = 10000
x = lapply(FUN = function(x) {y = cumsum(sample(c(-10, 10),size = x, TRUE));rands = abs(rnorm(n = x,mean = 0,sd = 1));return(list(High= y + abs(min(y))+rands,Low =  y + abs(min(y))-rands, vol.high = runif(1,0,1000),vol.low = runif(1,0,1000) ))},X = rep(list(x = paramsim$n),paramsim$n.curs))
results = array(0,c(N,2222))
for(i in 1:N)
{
  #paramsim = list(seed = i*10, n = 10000, n.curs = 500,buypercent=0.9)  
  #parambot$balance=balance
  paramsim = list(seed = i*10, n = 10, n.curs = 50,buypercent=1)
  parambot = list(limit_trade_pairs  =runif(1,1,paramsim$n.curs), min_buy_size = runif(1,0,100), ban_days_new_pair = runif(1,0,paramsim$n),days_to_keep = runif(1,0,paramsim$n), ban_days = runif(1,0,paramsim$n), min_profit = runif(1,0.0001,1), sell_loss = runif(1,-1,-0.0001), min_currency_price = runif(1,0.000001,100),  trade_step =  runif(1,1,paramsim$n-1), balance = balance, x.hist = x)
  res = simulateAction(paramsim,parambot)
  results[i,1111] = res$balance -  balance 
  balance = res$balance
  #compute some technical indicators from the previous period to reduce dimensionality of the data
  
    #adxs = unlist(lapply(FUN = function(X) { TTR::ADX(HLC = cbind(X$High,X$Low,(X$High+X$Low)/2))[paramsim$n]},X = x))
    #aroons = unlist(lapply(FUN = function(X) { TTR::aroon(HL = cbind(X$High,X$Low))[paramsim$n]},X = x))
    #atrs = unlist(lapply(FUN = function(X) { TTR::ATR(HLC = cbind(X$High,X$Low,(X$High+X$Low)/2))[paramsim$n]},X = x))
    #bbands = unlist(lapply(FUN = function(X) { TTR::BBands(HLC = cbind(X$High,X$Low,(X$High+X$Low)/2))[paramsim$n]},X = x))
    #ccis = unlist(lapply(FUN = function(X) { TTR::CCI(HLC = cbind(X$High,X$Low,(X$High+X$Low)/2))[paramsim$n]},X = x))
    #chaikinADS = unlist(lapply(FUN = function(X) { TTR::chaikinAD(HLC = cbind(X$High,X$Low,(X$High+X$Low)/2),volume = x[[1]]$vol.high)[paramsim$n]},X = x))
    
  #add chain's history to the arrays
  results[i,10:1109] = unlist(x)
  x = res$x
  results[i,1] = parambot$limit_trade_pairs
  results[i,2] = parambot$min_buy_size
  results[i,3] = parambot$ban_days_new_pair
  results[i,4] = parambot$days_to_keep
  results[i,5] = parambot$ban_days
  results[i,6] = parambot$min_profit
  results[i,7] = parambot$sell_loss
  results[i,8] = parambot$min_currency_price
  results[i,9] = parambot$trade_step
  results[i,1110] = res$balance
  
  if(i>1)
    results[i,1112:2222]=results[i-1,1:1111]
  print(results[i,1111])
  
  #plot chains' history if required
  #plot(x[[1]]$High,col=2)
  #points(x[[1]]$Low,col=3)
  
}
results=results[-1,]
N=N-1
#random bot is a stupid bot
plot(results[,1110])
plot(results[,1111])

#naively select the set of hyperparameters maximizing the reward
i = which(results[,1111]==max(results[,1111]))
 parambot = list(limit_trade_pairs  =runif(1,1,paramsim$n.curs), min_buy_size = runif(1,0,100), ban_days_new_pair = runif(1,0,paramsim$n),days_to_keep = runif(1,0,paramsim$n), ban_days = runif(1,0,paramsim$n), min_profit = runif(1,0.0001,1), sell_loss = runif(1,-1,-0.0001), min_currency_price = runif(1,0.000001,100),  trade_step =  runif(1,1,paramsim$n-1), balance = balance, x.hist = x)
 parambot$limit_trade_pairs  = results[i,1]
 parambot$min_buy_size =  results[i,2]
 parambot$ban_days_new_pair =  results[i,3]
 parambot$days_to_keep =  results[i,4]
 parambot$ban_days  =  results[i,5]
 parambot$min_profit  =  results[i,6] 
 parambot$sell_loss  =  results[i,7]
 parambot$min_currency_price  =  results[i,8] 
 parambot$trade_step  =  results[i,9]
 
 #see what happens with the bot
 N=15001
 balance = 10000
 x = lapply(FUN = function(x) {y = cumsum(sample(c(-10, 10),size = x, TRUE));rands = abs(rnorm(n = x,mean = 0,sd = 1));return(list(High= y + abs(min(y))+rands,Low =  y + abs(min(y))-rands, vol.high = runif(1,0,1000),vol.low = runif(1,0,1000) ))},X = rep(list(x = paramsim$n),paramsim$n.curs))
 results1 = array(0,c(N,2222))
 for(i in 1:N)
 {
   #paramsim = list(seed = i*10, n = 10000, n.curs = 500,buypercent=0.9)  
   #parambot$balance=balance
   paramsim = list(seed = i*10, n = 10, n.curs = 50,buypercent=1)
   
   res = simulateAction(paramsim,parambot)
   results1[i,1111] = res$balance -  balance 
   balance = res$balance
   #compute some technical indicators from the previous period to reduce dimensionality of the data
   parambot$x = x
   #adxs = unlist(lapply(FUN = function(X) { TTR::ADX(HLC = cbind(X$High,X$Low,(X$High+X$Low)/2))[paramsim$n]},X = x))
   #aroons = unlist(lapply(FUN = function(X) { TTR::aroon(HL = cbind(X$High,X$Low))[paramsim$n]},X = x))
   #atrs = unlist(lapply(FUN = function(X) { TTR::ATR(HLC = cbind(X$High,X$Low,(X$High+X$Low)/2))[paramsim$n]},X = x))
   #bbands = unlist(lapply(FUN = function(X) { TTR::BBands(HLC = cbind(X$High,X$Low,(X$High+X$Low)/2))[paramsim$n]},X = x))
   #ccis = unlist(lapply(FUN = function(X) { TTR::CCI(HLC = cbind(X$High,X$Low,(X$High+X$Low)/2))[paramsim$n]},X = x))
   #chaikinADS = unlist(lapply(FUN = function(X) { TTR::chaikinAD(HLC = cbind(X$High,X$Low,(X$High+X$Low)/2),volume = x[[1]]$vol.high)[paramsim$n]},X = x))
   
   #add chain's history to the arrays
   results1[i,10:1109] = unlist(x)
   x = res$x
   results1[i,1] = parambot$limit_trade_pairs
   results1[i,2] = parambot$min_buy_size
   results1[i,3] = parambot$ban_days_new_pair
   results1[i,4] = parambot$days_to_keep
   results1[i,5] = parambot$ban_days
   results1[i,6] = parambot$min_profit
   results1[i,7] = parambot$sell_loss
   results1[i,8] = parambot$min_currency_price
   results1[i,9] = parambot$trade_step
   results1[i,1110] = res$balance
   
   if(i>1)
     results1[i,1112:2222]=results1[i-1,1:1111]
   print(results1[i,1111])
   
   #plot chains' history if required
   #plot(x[[1]]$High,col=2)
   #points(x[[1]]$Low,col=3)
   
 }
 results1=results1[-1,]
 N=N-1

#naively trained bot is quite a bit better
plot(results1[,1110])
plot(results1[,1111]) 

#final profit, well this one is much smarter than the random bot, earns some money for us even
results1[N,1110]-results1[1,1110]-10000
#see the profit at all time points

plot(results1[,1110]-results1[1,1110]-10000)
lines(y=rep(mean(results1[,1110]-results1[1,1110]-10000),N),x = (1:N),col=2)
lines(y=rep(0,N),x = (1:N),col=3)

#cool, but now make a step further and try the Bayesian bot suggested by Hubin et. all 2017 out



#select all of the hyperparameters making the positive reward in the initial simulations
eps = 50000
ids = which(results[,1111]>eps)
lids = length(ids)
#define the uniform prior probability for them to be addressed
probs = array(1,lids)

#see what happens with the Bayesian bot
N=15001
balance = 10000
x = lapply(FUN = function(x) {y = cumsum(sample(c(-10, 10),size = x, TRUE));rands = abs(rnorm(n = x,mean = 0,sd = 1));return(list(High= y + abs(min(y))+rands,Low =  y + abs(min(y))-rands, vol.high = runif(1,0,1000),vol.low = runif(1,0,1000) ))},X = rep(list(x = paramsim$n),paramsim$n.curs))
results1 = array(0,c(N,2222))
for(i in 1:N)
{
  #paramsim = list(seed = i*10, n = 10000, n.curs = 500,buypercent=0.9)  
  #parambot$balance=balance
  
  paramsim = list(seed = i*10, n = 10, n.curs = 50,buypercent=1)
  parambot$x = x
  
  #now choose the posterior distribution's paraamteters of the grid
  ii = sample.int(n = lids,size = 1,prob = probs)
  parambot$limit_trade_pairs  = results[ii,1]
  parambot$min_buy_size =  results[ii,2]
  parambot$ban_days_new_pair =  results[ii,3]
  parambot$days_to_keep =  results[ii,4]
  parambot$ban_days  =  results[ii,5]
  parambot$min_profit  =  results[ii,6] 
  parambot$sell_loss  =  results[ii,7]
  parambot$min_currency_price  =  results[ii,8] 
  parambot$trade_step  =  results[ii,9]
  
  res = simulateAction(paramsim,parambot)
  results1[i,1111] = res$balance -  balance
  #update the posterior of the bot's distribution of hyperparameters
  if(results1[i,1111]>0)
  {
    probs[ii]=probs[ii]+1
  }else if(results1[i,1111]<0)
  {
    probs[ii]=probs[ii]-1
    if(probs[ii]==0)
      probs[ii]=1
  }
  balance = res$balance
  
  #compute some technical indicators from the previous period to reduce dimensionality of the data
  
  #adxs = unlist(lapply(FUN = function(X) { TTR::ADX(HLC = cbind(X$High,X$Low,(X$High+X$Low)/2))[paramsim$n]},X = x))
  #aroons = unlist(lapply(FUN = function(X) { TTR::aroon(HL = cbind(X$High,X$Low))[paramsim$n]},X = x))
  #atrs = unlist(lapply(FUN = function(X) { TTR::ATR(HLC = cbind(X$High,X$Low,(X$High+X$Low)/2))[paramsim$n]},X = x))
  #bbands = unlist(lapply(FUN = function(X) { TTR::BBands(HLC = cbind(X$High,X$Low,(X$High+X$Low)/2))[paramsim$n]},X = x))
  #ccis = unlist(lapply(FUN = function(X) { TTR::CCI(HLC = cbind(X$High,X$Low,(X$High+X$Low)/2))[paramsim$n]},X = x))
  #chaikinADS = unlist(lapply(FUN = function(X) { TTR::chaikinAD(HLC = cbind(X$High,X$Low,(X$High+X$Low)/2),volume = x[[1]]$vol.high)[paramsim$n]},X = x))
  
  #add chain's history to the arrays
  results1[i,10:1109] = unlist(x)
  x = res$x
  results1[i,1] = parambot$limit_trade_pairs
  results1[i,2] = parambot$min_buy_size
  results1[i,3] = parambot$ban_days_new_pair
  results1[i,4] = parambot$days_to_keep
  results1[i,5] = parambot$ban_days
  results1[i,6] = parambot$min_profit
  results1[i,7] = parambot$sell_loss
  results1[i,8] = parambot$min_currency_price
  results1[i,9] = parambot$trade_step
  results1[i,1110] = res$balance
  
  if(i>1)
    results1[i,1112:2222]=results1[i-1,1:1111]
  print(results1[i,1111])
  
  #plot chains' history if required
  #plot(x[[1]]$High,col=2)
  #points(x[[1]]$Low,col=3)
  
}
results1=results1[-1,]
N=N-1

#Bayesian bot is even smarter but yet makes mistakes now and then. Oh, boy...
plot(results1[,1110])
plot(results1[,1111]) 

#final profit
results1[N,1110]-results1[1,1110]-10000

#see the profit at all time points

plot(results1[,1110]-results1[1,1110]-10000)
lines(y=rep(mean(results1[,1110]-results1[1,1110]-10000),N),x = (1:N),col=2)
lines(y=rep(0,N),x = (1:N),col=3)

#!!!but let us try the smart Deeply Trained bot out!!! the rest is not implemented yet!?

#!!!!!!!!!!the code below is not working yet!!!!!!!!!!!


#a linear stack of layers
model<-keras_model_sequential()
#configuring the Model
# model %>%
#   layer_conv_1d(filter=100 ,kernel_size=c(5),input_shape = c(500,1))%>%  
#   layer_activation("relu") %>%
#   #dropout layer to avoid overfitting
#   layer_dense(512) %>%
#   layer_max_pooling_1d(pool_size = c(3)) %>%
#   layer_dropout(0.1) %>%
#   layer_flatten() %>%  
#   layer_dense(512) %>%
#   layer_activation("relu") %>%  
#   layer_dropout(0.5) %>%  
#   #output layer-10 classes-10 units  
#   layer_dense(1) %>%
#   layer_activation("linear") 

model %>% 
  layer_dense(units = 2221, input_shape = 2221) %>% 
  layer_activation(activation = 'sigmoid') %>% 
  layer_dense(units = 1000) %>% 
  layer_activation(activation = 'sigmoid')%>% 
  layer_dropout(rate=0.4)%>%
  layer_dense(units = 100) %>% 
  layer_activation(activation = 'sigmoid')%>% 
  layer_dense(units = 1000) %>% 
  layer_activation(activation = 'sigmoid')%>% 
  layer_dense(units = 1) %>% 
  layer_activation("linear") 

#Model's Optimizer
#defining the type of optimizer-ADAM-Adaptive Momentum Estimation
opt<-optimizer_adam( lr= 0.0001 , decay = 1e-6 )
#lr-learning rate , decay - learning rate decay over each update

model %>%
compile(loss="mae",
optimizer=opt,metrics = "mae")
#Summary of the Model and its Architecture
summary(model)
ids = sample(x = N,size = (N)/2,replace = F)
train_x=results[ids,-1111]
train_y=results[ids,1111]
test_x=results[-ids,-1111]
test_y=results[-ids,1111]
  
#TRAINING PROCESS OF THE MODEL
model %>% fit( train_x,train_y,
                 epochs=100,batch_size = 1000, validation_data = list(test_x, test_y),
                 shuffle=TRUE)

MSE.test = model %>% evaluate(test_x, test_y)
MSE.train = model %>% evaluate(train_x,train_y)

#mean prediction 
my = mean(train_y)
mae = sum(abs(train_y - my))/length(train_y)


# make an Q-learning based bot
source("Agent.R")
source("Memory.R")
source("train.R")

N_EPISODE <- 500
BATCH_SIZE <- 32

library(gym)


tf$reset_default_graph()

agent <-
  Agent$new(
    input_shape = 9,
    output_dim = 2,
    epsilon_last_episode = 100
  )

memory <- Memory$new(capacity = 50000)

rewards <- c()

with(tf$Session() %as% sess, {
  init <- tf$global_variables_initializer()
  sess$run(init)
  #start Q learning 
  for (episode_i in 1:N_EPISODE) {
    done <- FALSE
    #reinitialize the stock markets
    s <- env$reset()
    total_reward = 0
    
    while (!done) {
      #decide on the set of bot's paramaters given the market sition history and paramaters' history
      a <- agent$get_action(state_ = s, step = episode_i)
      #get the reward from action a
      ret <- env$step(action = a)
      
      s2 <- ret[[1]]#the new state
      r <- ret[[2]]#the current reward
      done <- ret[[3]]#if the node is terminal
      
      memory$push(s, a, r, done, s2)#update the POMDP history
      
      #update the agent via training
      if (memory$length > BATCH_SIZE) {
        batch <- memory$sample(BATCH_SIZE)
        train(agent, batch)
      }
      s <- s2#change the state
      total_reward <- total_reward + r#calculate total reward
    }
    
    cat(
      sprintf(
        "[Episode: %4d] Reward: %4d, Epsilon: %.3f\n",
        episode_i,
        total_reward,
        agent$epsilon
      )
    )
    
    rewards <- append(rewards, total_reward)
    
    if (length(rewards) > 100) {
      rewards <- rewards[2:length(rewards)]
      
      if (mean(rewards) > 195) {
        cat("Game Cleared")
        break
      }
    }
  }
})

env$close()
  

policy = EpsGreedyQPolicy()
memory = SequentialMemory(limit=50000, window_length=1)
dqn = DQNAgent(model=model, nb_actions=nb_actions, memory=memory, nb_steps_warmup=10,
               target_model_update=1e-2, policy=policy)
dqn.compile(Adam(lr=1e-3), metrics='mae')
            
            


env <- gridworldEnvironment
print(env)

# Define state and action sets
states <- c("s1", "s2", "s3", "s4")
actions <- c("up", "down", "left", "right")

# Sample N = 1000 random sequences from the environment
data <- sampleExperience(N = 1000, env = env, states = states, actions = actions)
head(data)

# Define reinforcement learning parameters
control <- list(alpha = 0.1, gamma = 0.5, epsilon = 0.1)

# Perform reinforcement learning
model <- ReinforcementLearning(data, s = "State", a = "Action", r = "Reward", 
                               s_new = "NextState", control = control)

# Print result
print(model)

      
# Define reinforcement learning parameters
control <- list(alpha = 0.1, gamma = 0.5, epsilon = 0.1)

# Sample N = 1000 sequences from the environment using epsilon-greedy action selection
data_new <- sampleExperience(N = 1000, env = env, states = states, actions = actions, 
                             model = model, actionSelection = "epsilon-greedy", 
                             control = control)
head(data_new)

# Update the existing policy using new training data
model_new <- ReinforcementLearning(data_new, s = "State", a = "Action", r = "Reward", 
                                   s_new = "NextState", control = control, model = model)

# Print result
print(model_new)
