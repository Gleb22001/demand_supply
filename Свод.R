library(readxl)



####Продовльственные товары####
data1 <- read_xlsx("Продовольственные товары_декомпозиции.xlsx", sheet="Sheet1", col_names=FALSE)
data2 <- read_xlsx("Продовольственные товары_декомпозиции.xlsx", sheet="Sheet2", col_names=FALSE)
food_ipc <- read_xlsx("food_ipc.xlsx", sheet="ln")

q_demand <- data1[ ,1]
q_supply <- data1[ ,2]

p_demand <- data2[ ,1]
p_supply <- data2[ ,2]

for (i in 2:17) {
  d1=paste0("Sheet", 2*i-1)
  d2=paste0("Sheet", 2*i)
  data1 <- read_xlsx("Продовольственные товары_декомпозиции.xlsx", sheet=d1, col_names=FALSE)
  data2 <- read_xlsx("Продовольственные товары_декомпозиции.xlsx", sheet=d2, col_names=FALSE)
  q_demand <- cbind(q_demand, data1[ ,1])
  q_supply <- cbind(q_supply, data1[ ,2])
  p_demand <- cbind(p_demand, data2[ ,1])
  p_supply <- cbind(p_supply, data2[ ,2])
  
  
}

names(q_demand) <- names(food_ipc)[-1]
names(q_supply) <- names(food_ipc)[-1]


names(p_demand) <- names(food_ipc)[-1]
names(p_supply) <- names(food_ipc)[-1]

q_demand_prod <- q_demand
q_supply_prod <- q_supply

p_demand_prod <- p_demand
p_supply_prod <- p_supply

q_demand_prod$date <- food_ipc$...1[8:154]
q_supply_prod$date <- food_ipc$...1[8:154]


p_demand_prod$date <- food_ipc$...1[8:154]
p_supply_prod$date <- food_ipc$...1[8:154]

writexl::write_xlsx(q_demand, "q_demand_prod.xlsx")
writexl::write_xlsx(q_supply, "q_supply_prod.xlsx")

writexl::write_xlsx(p_demand, "p_demand_prod.xlsx")
writexl::write_xlsx(p_supply, "p_supply_prod.xlsx")



####Непроды####
data1 <- read_xlsx("Неродовольственные товары_декомпозиции.xlsx", sheet="Sheet1", col_names=FALSE)
data2 <- read_xlsx("Неродовольственные товары_декомпозиции.xlsx", sheet="Sheet2", col_names=FALSE)
sales_basic_nfood <- read_xlsx("sales_basic_nfood.xlsx", sheet="ln")


q_demand <- data1[ ,1]
q_supply <- data1[ ,2]


p_demand <- data2[ ,1]
p_supply <- data2[ ,2]


for (i in 2:12) {
  d1=paste0("Sheet", 2*i-1)
  d2=paste0("Sheet", 2*i)
  data1 <- read_xlsx("Неродовольственные товары_декомпозиции.xlsx", sheet=d1, col_names=FALSE)
  data2 <- read_xlsx("Неродовольственные товары_декомпозиции.xlsx", sheet=d2, col_names=FALSE)
  q_demand <- cbind(q_demand, data1[ ,1])
  q_supply <- cbind(q_supply, data1[ ,2])
  p_demand <- cbind(p_demand, data2[ ,1])
  p_supply <- cbind(p_supply, data2[ ,2])
}


names(q_demand) <- names(sales_basic_nfood)[-1]
names(q_supply) <- names(sales_basic_nfood)[-1]

names(p_demand) <- names(sales_basic_nfood)[-1]
names(p_supply) <- names(sales_basic_nfood)[-1]





q_demand_nprod <- q_demand
q_supply_nprod <- q_supply


p_demand_nprod <- p_demand
p_supply_nprod <- p_supply


q_demand_nprod$date <- food_ipc$...1[8:154]
q_supply_nprod$date <- food_ipc$...1[8:154]


p_demand_nprod$date <- food_ipc$...1[8:154]
p_supply_nprod$date <- food_ipc$...1[8:154]










writexl::write_xlsx(q_demand, "q_demand_nprod.xlsx")
writexl::write_xlsx(q_supply, "q_supply_nprod.xlsx")

writexl::write_xlsx(p_demand, "p_demand_nprod.xlsx")
writexl::write_xlsx(p_supply, "p_supply_nprod.xlsx")





####Услуги ####


data1 <- read_xlsx("Услуги_декомпозиции.xlsx", sheet="Sheet1", col_names=FALSE)
data2 <- read_xlsx("Услуги_декомпозиции.xlsx", sheet="Sheet2", col_names=FALSE)
services_ipc <- read_xlsx("services_ipc.xlsx", sheet="ln")


q_demand <- data1[ ,1]
q_supply <- data1[ ,2]


p_demand <- data2[ ,1]
p_supply <- data2[ ,2]


for (i in 2:13) {
  d1=paste0("Sheet", 2*i-1)
  d2=paste0("Sheet", 2*i)
  data1 <- read_xlsx("Услуги_декомпозиции.xlsx", sheet=d1, col_names=FALSE)
  data2 <- read_xlsx("Услуги_декомпозиции.xlsx", sheet=d2, col_names=FALSE)
  q_demand <- cbind(q_demand, data1[ ,1])
  q_supply <- cbind(q_supply, data1[ ,2])
  p_demand <- cbind(p_demand, data2[ ,1])
  p_supply <- cbind(p_supply, data2[ ,2])
  
}


names(q_demand) <- names(services_ipc)[-1]
names(q_supply) <- names(services_ipc)[-1]

names(p_demand) <- names(services_ipc)[-1]
names(p_supply) <- names(services_ipc)[-1]


q_demand_serv <- q_demand
q_supply_serv <- q_supply

p_demand_serv <- p_demand
p_supply_serv <- p_supply


q_demand_serv$date <- food_ipc$...1[8:154]
q_supply_serv$date <- food_ipc$...1[8:154]


p_demand_serv$date <- food_ipc$...1[8:154]
p_supply_serv$date <- food_ipc$...1[8:154]

writexl::write_xlsx(q_demand, "q_demand_serv.xlsx")
writexl::write_xlsx(q_supply, "q_supply_serv.xlsx")
writexl::write_xlsx(p_demand, "p_demand_serv.xlsx")
writexl::write_xlsx(p_supply, "p_supply_serv.xlsx")

####Регионы####

data1 <- read_xlsx("Регионы_декомпозиции.xlsx", sheet="Sheet1", col_names=FALSE)
data2 <- read_xlsx("Регионы_декомпозиции.xlsx", sheet="Sheet2", col_names=FALSE)

q_demand <- data1[ ,1]
q_supply <- data1[ ,2]


p_demand <- data2[ ,1]
p_supply <- data2[ ,2]


for (i in 2:8) {
  d1=paste0("Sheet", 2*i-1)
  d2=paste0("Sheet", 2*i)
  data1 <- read_xlsx("Регионы_декомпозиции.xlsx", sheet=d1, col_names=FALSE)
  data2 <- read_xlsx("Регионы_декомпозиции.xlsx", sheet=d2, col_names=FALSE)
  q_demand <- cbind(q_demand, data1[ ,1])
  q_supply <- cbind(q_supply, data1[ ,2])
  p_demand <- cbind(p_demand, data2[ ,1])
  p_supply <- cbind(p_supply, data2[ ,2])
}

names <- c("Москва", "Санкт-Петербург", "Краснодарский край", "Ставропольский край", 
           "Татарстан", "Тюменская область", "Красносярский край", "Саха")

names(q_demand) <- names
names(q_supply) <- names


names(p_demand) <- names
names(p_supply) <- names



writexl::write_xlsx(q_demand, "q_demand_region.xlsx")
writexl::write_xlsx(q_supply, "q_supply_region.xlsx")
writexl::write_xlsx(p_demand, "p_demand_region.xlsx")
writexl::write_xlsx(p_supply, "p_supply_region.xlsx")

