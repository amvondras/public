# Assess methylation bias along reads using bismark output

_Split files_

Each M-bias file need to be split into 6 files of equal length, with extra information omitted, so it can be imported to R for plotting
        
        head -n 156 CF01*.M-bias.txt | tail -n 154 | sed 's/ /_/g' | column -t | head

        #File order:
        #0. CpG context (R1)
        #1. CHG context (R1)
        #2. CHH context (R1)
        #3. CpG context (R2)
        #4. CHG context (R2)
        #5. CHH context (R2)

        for i in $( ls *.txt | grep -v "sub" | sed 's/.deduplicated.M-bias.txt//g' ); do  cat ${i}.deduplicated.M-bias.txt | grep -v "=" | grep -v "(" | sed 's/ /_/g' | awk -v env_var="$i" -v RS= '{print > (env_var "." NR-1 ".M-bias.sub.txt")}' ; done
        
_Import to R and plot_

      setwd("~/Desktop/Machine/methylome.project/methylation/bismark/on_cs/mbias/")
      
      data=list(
        CF = list(
          CF01 = list(cpg.r1 = read.delim("CF01.0.M-bias.sub.txt"), 
                       chg.r1 = read.delim("CF01.1.M-bias.sub.txt"),
                       chh.r1 = read.delim("CF01.2.M-bias.sub.txt"),
                       cpg.r2 = read.delim("CF01.3.M-bias.sub.txt"), 
                       chg.r2 = read.delim("CF01.4.M-bias.sub.txt"),
                       chh.r2 = read.delim("CF01.5.M-bias.sub.txt")),
          CF03 = list(cpg.r1 = read.delim("CF03.0.M-bias.sub.txt"), 
                      chg.r1 = read.delim("CF03.1.M-bias.sub.txt"),
                      chh.r1 = read.delim("CF03.2.M-bias.sub.txt"),
                      cpg.r2 = read.delim("CF03.3.M-bias.sub.txt"), 
                      chg.r2 = read.delim("CF03.4.M-bias.sub.txt"),
                      chh.r2 = read.delim("CF03.5.M-bias.sub.txt")),
          CF04 = list(cpg.r1 = read.delim("CF04.0.M-bias.sub.txt"), 
                      chg.r1 = read.delim("CF04.1.M-bias.sub.txt"),
                      chh.r1 = read.delim("CF04.2.M-bias.sub.txt"),
                      cpg.r2 = read.delim("CF04.3.M-bias.sub.txt"), 
                      chg.r2 = read.delim("CF04.4.M-bias.sub.txt"),
                      chh.r2 = read.delim("CF04.5.M-bias.sub.txt"))
        ),
        CS = list(
          CS06 = list(cpg.r1 = read.delim("CS06.0.M-bias.sub.txt"), 
                      chg.r1 = read.delim("CS06.1.M-bias.sub.txt"),
                      chh.r1 = read.delim("CS06.2.M-bias.sub.txt"),
                      cpg.r2 = read.delim("CS06.3.M-bias.sub.txt"), 
                      chg.r2 = read.delim("CS06.4.M-bias.sub.txt"),
                      chh.r2 = read.delim("CS06.5.M-bias.sub.txt")),
          CS08 = list(cpg.r1 = read.delim("CS08.0.M-bias.sub.txt"), 
                      chg.r1 = read.delim("CS08.1.M-bias.sub.txt"),
                      chh.r1 = read.delim("CS08.2.M-bias.sub.txt"),
                      cpg.r2 = read.delim("CS08.3.M-bias.sub.txt"), 
                      chg.r2 = read.delim("CS08.4.M-bias.sub.txt"),
                      chh.r2 = read.delim("CS08.5.M-bias.sub.txt")),
          CS47 = list(cpg.r1 = read.delim("CS47.0.M-bias.sub.txt"), 
                      chg.r1 = read.delim("CS47.1.M-bias.sub.txt"),
                      chh.r1 = read.delim("CS47.2.M-bias.sub.txt"),
                      cpg.r2 = read.delim("CS47.3.M-bias.sub.txt"), 
                      chg.r2 = read.delim("CS47.4.M-bias.sub.txt"),
                      chh.r2 = read.delim("CS47.5.M-bias.sub.txt"))
        ),
        SB = list(
          SB01 = list(cpg.r1 = read.delim("SB01.0.M-bias.sub.txt"), 
                      chg.r1 = read.delim("SB01.1.M-bias.sub.txt"),
                      chh.r1 = read.delim("SB01.2.M-bias.sub.txt"),
                      cpg.r2 = read.delim("SB01.3.M-bias.sub.txt"), 
                      chg.r2 = read.delim("SB01.4.M-bias.sub.txt"),
                      chh.r2 = read.delim("SB01.5.M-bias.sub.txt")),
          SB06 = list(cpg.r1 = read.delim("SB06.0.M-bias.sub.txt"), 
                      chg.r1 = read.delim("SB06.1.M-bias.sub.txt"),
                      chh.r1 = read.delim("SB06.2.M-bias.sub.txt"),
                      cpg.r2 = read.delim("SB06.3.M-bias.sub.txt"), 
                      chg.r2 = read.delim("SB06.4.M-bias.sub.txt"),
                      chh.r2 = read.delim("SB06.5.M-bias.sub.txt")),
          SB14 = list(cpg.r1 = read.delim("SB14.0.M-bias.sub.txt"), 
                      chg.r1 = read.delim("SB14.1.M-bias.sub.txt"),
                      chh.r1 = read.delim("SB14.2.M-bias.sub.txt"),
                      cpg.r2 = read.delim("SB14.3.M-bias.sub.txt"), 
                      chg.r2 = read.delim("SB14.4.M-bias.sub.txt"),
                      chh.r2 = read.delim("SB14.5.M-bias.sub.txt"))
        )
      )  

        for (l1 in 1:length(data)){
          print(names(data)[l1])
          for(l2 in 1:length(data[[l1]])){
            print(names(data[[l1]])[l2])
            for(l3 in 1:length(data[[l1]][[l2]])){
              print(names(data[[l1]][[l2]])[l3])
              data[[l1]][[l2]][[l3]]$sample <- names(data[[l1]])[l2]
              data[[l1]][[l2]][[l3]]$set <- gsub(".r.", "", names(data[[l1]][[l2]])[l3])
              data[[l1]][[l2]][[l3]]$read <- gsub(".*\\.r", "R", names(data[[l1]][[l2]])[l3])
            }
          }
        }; head(data$CF$CF01$cpg.r1)

        # Collapse list into single dataframe
              my.df <- data.frame()  
              for (l1 in 1:length(data)){
                print(names(data)[l1])
                for(l2 in 1:length(data[[l1]])){
                  print(names(data[[l1]])[l2])
                  for(l3 in 1:length(data[[l1]][[l2]])){
                    print(head(data[[l1]][[l2]][[l3]]))
                    my.df <- rbind(my.df, data.frame(data[[l1]][[l2]][[l3]]))
                  }
                }
              }; rm(l1, l2, l3); head(my.df)

        table(my.df$sample, my.df$set, my.df$read)  

        ggplot(data=my.df, aes(x=position, y=X._methylation, color=set)) + 
          geom_point(size=0.5) + 
          xlab("Read position") + 
          ylab("% Methylation") +
          theme_classic() +
          facet_grid(sample~read) +
          scale_x_continuous(breaks=c(0,7, 25, 50,75,100,150,155)) +
          geom_vline(xintercept = 7, size=0.1) + 
          geom_vline(xintercept = 150, size=0.1) +  
          geom_vline(xintercept = 100, size=0.1) 

It looks like the first 15 bases and last two bases of R1 and R2 should be ignored by bismark. 
        
[mbias.pdf](https://github.com/amvondras/trio/files/6594707/mbias.pdf)

<img width="559" alt="Screen Shot 2021-06-03 at 5 28 17 PM" src="https://user-images.githubusercontent.com/33852065/120728159-241b4100-c491-11eb-882c-df95a12d749f.png">

After removing the first 7 bases of R1 and R2 and the last two bases of R2:

![image](https://user-images.githubusercontent.com/33852065/120737991-ed025b00-c4a3-11eb-957d-89c737d9d7e1.png)

After removing the first 15 bases and last two bases of R1 and R2:
![image](https://user-images.githubusercontent.com/33852065/121131545-cf572d80-c7e4-11eb-9a31-57bf6637a67b.png)



