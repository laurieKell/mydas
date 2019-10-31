stockDef=rbind(expand.grid(stock="codns",division=c("4A","4B","4C", "7D", "3A", "4"), stock="codns"),
               expand.grid(stock="hadns",division=c("4A","4B","4C","3A","4"), "nseahad"),        
               expand.grid(stock="whgns",division=c("4A","4B","4C","4","7D"), "nseawhg"),        
               expand.grid(stock="pokns",division=c("4A","4B","4C","4","6A","6B","3A"), "nseapok"),        
               expand.grid(stock="plens",division=c("4A","4B","4C","4", "3A"), "nseaple"),        
               expand.grid(stock="solns",division=c("4A","4B","4C","4"), "nseasol"),
               expand.grid(stock="codcs",division=c("7E", "7F","7H", "7G","7J", "7K"), "cseacod"),
               expand.grid(stock="hadcs",division=c("7B","7C", "7D","7E", "7F","7H", "7G","7J","7K"), "cseahad"),
               expand.grid(stock="whgcs",division=c("7B","7C","7H", "7E", "7F","7H", "7G","7J","7K")))