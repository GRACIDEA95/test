##############実験用
aList <- JP.ID.aList
pList <- JP.ID.pList

aList <- JP.GF.aList
pList <- JP.GF.pList

aList <- list()
for (i in 1:16){
  aList[[i]] <- JP.GF.aList[[i]]
}
for (i in 1:4){
  aList[[16+i]] <- JP.GF.aList[[17+i]]
}
pList <- list()
for (i in 1:16){
  pList[[i]] <- JP.GF.pList[[i]]
}
for (i in 1:4){
  pList[[16+i]] <- JP.GF.pList[[17+i]]
}

a <- tkgetOpenFile(filetypes = "{{CSV Files} {.csv}}")
  tclvalue(a)
  mutation_df <- read.csv(tclvalue(a))
  mutation_paternal <- mutation_df[,2:6]
  mutation_maternal <- mutation_df[,8:12]

pedigree <- rbind(c(0,0),c(0,0),c(1,2),c(0,0),c(1,2),c(3,4),c(3,4)) #家系図行列
pedigree <- rbind(c(0,0),c(0,0),c(0,0),c(1,2),c(1,2),c(2,3)) #家系図行列
sample_one <- Pedigree.make_linkage_mutation(pedigree,1,aList,pList,mutation_paternal,mutation_maternal,linkage.loci,RR)

t <- proc.time()
LR <- PES_function_linkage(names,pedigree,sample_one,target_number,0,mutation_paternal,mutation_maternal,aList,pList,linkage.loci,RR)
proc.time() - t

t <- proc.time()
LR <- PES_function(names,pedigree,sample_one,target_number,0,mutation_paternal,mutation_maternal,aList,pList)
proc.time() - t

####エラーサンプル20180807
error_sample <- sample_one 

  window <- tktoplevel()
  frame <- tkframe(window)

#######################ゴミ捨て場

  simulation_item <- tclVar("")
  simulation_mean <- tclVar("")
  simulation_median <- tclVar("")
  simulation_SD <- tclVar("")
  simulation_LR500 <- tclVar("")
  simulation_LR50 <- tclVar("")
  simulation_LR5 <- tclVar("")
s <- tclVar("disable")
s <- tclVar("normal")
bu <- tkbutton(tab.Files,text="button",state=tclvalue(s))
tkpack(bu)
tkconfigure(bu,state="normal")

  combobox_Linkage.locus1[[1]] <- ttkcombobox(frame.Linkage.select.input[[1]],values=var_locus,state="readonly",textvariable=locus1,font="our_font",width="10")
  combobox_Linkage.locus2[[1]] <- ttkcombobox(frame.Linkage.select.input[[1]],values=var_locus,state="readonly",textvariable=locus2,font="our_font",width="10")
  entry_Linkage.RR[[1]] <- tkentry(frame.Linkage.select.input[[1]], textvariable=RR, width="10", font="our_font", background="#ffffff")

  button_Linkage.add <- tkbutton(frame.Linkage.button,text="追加",font="our_font",width="10",
                                 command = function(locus1,locus2){
                                   linkage_number <<- linkage_number + 1

                                   frame.Linkage.select.input[[linkage_number]] <<- tkframe(frame.Linkage.select)
                                   combobox_Linkage.locus1[[linkage_number]] <- ttkcombobox(frame.Linkage.select.input[[linkage_number]],values=var_locus,state="readonly",textvariable=locus1[linkage_number],font="our_font",width="10")
                                   combobox_Linkage.locus2[[linkage_number]] <- ttkcombobox(frame.Linkage.select.input[[linkage_number]],values=var_locus,state="readonly",textvariable=locus2[linkage_number],font="our_font",width="10")
                                   entry_Linkage.RR[[linkage_number]] <- tkentry(frame.Linkage.select.input[[linkage_number]], textvariable=RR[linkage_number], width="10", font="our_font", background="#ffffff")

                                   tkpack(frame.Linkage.select.input[[linkage_number]],fill="none",anchor="w",padx="25",pady="10")
                                   tkpack(combobox_Linkage.locus1[[linkage_number]],padx="15",side="left")
                                   tkpack(combobox_Linkage.locus2[[linkage_number]],padx="15",side="left")
                                   tkpack(entry_Linkage.RR[[linkage_number]],padx="15",side="left")
                                                     }
                                 )
  button_Linkage.vanish <- tkbutton(frame.Linkage.button,text="1つ削除",font="our_font",width="10",
                                    command = function(){
                                      if (linkage_number > 1){
                                        tkdestroy(frame.Linkage.select.input[[linkage_number]])
                                        linkage_number <<- linkage_number - 1
                                      }
                                                        }
                                    )

  button_Linkage.next <- tkbutton(frame.Linkage.button,text="次へ",font="our_font",width="10",
                                    command = function(){
                                      LR_result <<- Calculate(filepath_LR.frequency1,filepath_LR.DNAtype1,filepath_LR.mutation1,as.integer(tclvalue(tkget(combobox_range1))))
                                      if (length(LR_result) > 1){
                                        Log(LR_result,log_LR_result)
                                        tkpack(label_LR_result)
                                        tk2notetab.select(Tabs, "LR Calculation")
                                      }
                                                        }
                                    )
  tkpack(frame.Linkage,fill="none")
  tkpack(frame.Linkage.label,fill="none",anchor="w",padx="25",pady="25")
  tkpack(frame.Linkage.select,fill="none")
  tkpack(frame.Linkage.select.input[[1]],fill="none",anchor="w",padx="25",pady="10")
  tkpack(frame.Linkage.button,fill="none",anchor="w",padx="25",pady="25")

  tkpack(label_Linkage.locus1,fill="none",side="left",padx="15")
  tkpack(label_Linkage.locus2,fill="none",side="left",padx="15")
  tkpack(label_Linkage.RR,fill="none",anchor="w",side="left",padx="15")
  tkpack(combobox_Linkage.locus1[[1]],padx="15",side="left")
  tkpack(combobox_Linkage.locus2[[1]],padx="15",side="left")
  tkpack(entry_Linkage.RR[[1]],padx="15",side="left")

  tkpack(button_Linkage.add,padx="15",side="left")
  tkpack(button_Linkage.vanish,padx="15",side="left")
  tkpack(button_Linkage.next,padx="15",side="left")




  tkinsert(entry_Linkage.RR[[1]],0,1)

        if (List[[1]][List[[2]][i,1]] != "x"){
          banish_number <- c(banish_number,List[[2]][i,1])
        }
        if (List[[1]][List[[2]][i,2]] != "x"){
          banish_number <- c(banish_number,List[[2]][i,2])
        }

      DNAtype.temp <- List[[3]][,,-b,]
      List[[2]][b:dim(List[[2]])[1],] <- apply(List[[2]][b:dim(List[[2]])[1],],c(1,2),function(x){if (x)return(names[which(names == x) - 1])})

  frame1_target <- tkframe(frame1)

  values_target.Hp <- c("0")
  values_target.Hd <- c("0")
  var_target.Hp <- tclVar(values_target.Hp[1])
  var_target.Hd <- tclVar(values_target.Hd[1])
  label_target.Hp <- tklabel(frame1_target,text="Hp:",font="our_font",width="4",anchor="w")
  combobox_target.Hp <- ttkcombobox(frame1_target, values=values_target.Hp, state="disable", textvariable=var_target.Hp, font="our_font",width="8")
  label_target.Hd <- tklabel(frame1_target,text="Hd:",font="our_font",width="4",anchor="w")
  combobox_target.Hd <- ttkcombobox(frame1_target, values=values_target.Hd, state="disable", textvariable=var_target.Hd, font="our_font",width="8")

  tkpack(frame1_target,fill="none",anchor="w",pady="1")


  tkpack(label_target.Hp,fill="none",anchor="w",side="left")
  tkpack(combobox_target.Hp,side="left")
  tkpack(label_target.Hd,fill="none",anchor="w",side="left")
  tkpack(combobox_target.Hd)

 tkconfigure(button_pedigree.check,state="disable")

  a <- "hello"
  b <- tclVar("")
  t <- tkentry(frame,textvariable=b,width="8",background="#ffffff")
  tkpack(frame)
  tkpack(t)
  tkentry(t)
  tclvalue(b)
tkText_widget.get

  frame2 <- tkframe(frame)
  frame2_1_1 <- tkframe(frame1)
  frame2_1_2 <- tkframe(frame1)
  frame2_2_1 <- tkframe(frame1)
  frame2_2_2 <- tkframe(frame1)
  frame2_3_1 <- tkframe(frame1)
  frame2_3_2 <- tkframe(frame1)
  frame2_4 <- tkframe(frame1)
  frame2_calculate <- tkframe(frame1)
  frame2_result <- tkframe(frame1)



  #log_result1 <- tclVar("ここに結果が表示されます")
      #tclvalue(log) <- "計算しています"


a <- tkgetOpenFile(filetypes = "{{CSV Files} {.csv}}")
  aplist <- read.csv(tclvalue(a))
  aplist[1,2]
  names(aplist)
  aplist[,2][!is.na(aplist[,2])]
  aplist[,1][!is.na(aplist[,2])]

  aList <- list()
  pList <- list()
  for (i in 1:(length(names(aplist))-1)){
    aList[[i]] <-  aplist[,1][!is.na(aplist[,i+1])]
    pList[[i]] <-  aplist[,i+1][!is.na(aplist[,i+1])]
  }

  var <- tclVar("hoge")
  label_LR.frequency <- tklabel(window,textvariable=var,font="our_font")
  tclvalue(var) <- "hogehoge"

  frame1 <- tkframe(window,bg="black")
  tclvalue(OK_LR.frequency) <- "○"

  tkdestroy(frame1)
  entry_range <- tkentry(frame1_4,text="hogehoge",font="our_font",width="4")
  
  tclvalue(tkget(combobox_menu))
  comboboxSelected()

  FileName <- tkgetOpenFile(filetypes = "{{CSV Files} {.csv}}")
  mutation_df <- read.csv(tclvalue(FileName))


  #msgbox_msg <- tkmessage(window,text="hogehogepiyppiyp",font="our_font",width=250)
  #tkpack(msgbox_msg)

  #tkconfigure(label_LR.frequency)
  #tkconfigure(button_LR.frequency,style="our_font")

  frequency_matrix <- read.csv("C:/Users/Hide/Documents/1_university/grade_4/R/PES/pList.csv")
  allelename <- names(frequency_matrix)
  paste(names(frequency_matrix))

  DNAtype_matrix <- read.csv("C:/Users/Hide/Documents/1_university/grade_4/R/PES/allele_sample.csv")
  as.character(DNAtype_matrix[,"name"][!is.na(DNAtype_matrix[,"name"])])
  cbind(DNAtype_matrix[,"paternal"][!is.na(DNAtype_matrix[,"paternal"])],DNAtype_matrix[,"maternal"][!is.na(DNAtype_matrix[,"maternal"])])
  DNAtype <- array(0,dim=c(2,length(names(frequency_matrix))-1,dim(DNAtype_matrix)[1]/2,1))
  DNAtype <- array(0,dim=c(2,15,7,1))
  for(i in 1:(dim(DNAtype_matrix)[1]/2)){
    DNAtype[1,,i,1] <- unlist(DNAtype_matrix[,-(1:3)][2 * i - 1,])
    DNAtype[2,,i,1] <- unlist(DNAtype_matrix[,-(1:3)][2 * i,])
  }
mode(DNAtype_matrix)
DNAtype_matrix[,-(1:3)][9,]
unlist(DNAtype_matrix[,-(1:3)][8,])


  window2 <- tktoplevel()
  tkfocus(window2)
  tkwm.title(window, "LRめいかー")

  a <- rbind(c(0,0),c(0,0),c(1,2),c(1,2),c(0,0),c(1,2),c(1,2),c(4,5),c(6,7))
  sample2 <- pedigree2(a,1,JP.ID.aList,JP.ID.pList)

  GetCursor()
cursor
window <- tktoplevel()
mes <- tklabel(window,text="hogehoge",font="our_font",cursor="xterm",width="250",justify="left",relief="groove",borderwidth="2")
tkpack(mes)

window <- tktoplevel()
mes <- ttkentry(window,font="our_font",width="250")
tkpack(mes)

a <- "11/22/3"
strsplit(a,"/")[[1]][length(strsplit(a,"/")[[1]])]

a <- c(1,4,2,5,4)
a2 <- tclVar("touhoku Line")

b <- c(2,3,4,5,6)
list(a,b)

write(a,file=Dir)
Dir <- "C:/Users/Hide/Documents/1_university/grade_4/R/GUI/test.txt"
b <- tkgetSaveFile(a2)

tclvalue(filepath_LR.frequency1) == "C:/Users/Hide/Documents/1_university/grade_4/R/PES/JP.GF.pList.csv"
tclvalue(filepath_LR.DNAtype1)
tclvalue(filepath_LR.mutation1)

fp_LR.frequency,fp_LR.DNAtype,fp_LR.mutation
filepath_LR.frequency1 <- tclVar("C:/Users/Hide/Documents/1_university/grade_4/R/PES/pList.csv")
filepath_LR.DNAtype1 <- tclVar("C:/Users/Hide/Documents/1_university/grade_4/R/PES/allele_sample.csv")
filepath_LR.mutation1 <- tclVar("C:/Users/Hide/Documents/1_university/grade_4/R/PES/mutation_rates.csv")
fp_LR.frequency <- tclVar("C:/Users/Hide/Documents/1_university/grade_4/R/PES/JP.GF.pList.csv")
fp_LR.DNAtype <- tclVar("C:/Users/Hide/Documents/1_university/grade_4/R/PES/allele_sample3_H1.csv")
fp_LR.DNAtype <- tclVar("C:/Users/Hide/Documents/1_university/grade_4/R/PES/allele_sample3_H1H2.csv")
fp_LR.DNAtype <- tclVar("C:/Users/Hide/Documents/1_university/grade_4/R/PES/allele_sample4_H1H2.csv")
fp_LR.DNAtype <- tclVar("C:/Users/Hide/Documents/1_university/grade_4/R/PES/allele_sample5_H1H2.csv")
fp_LR.DNAtype <- tclVar("C:/Users/Hide/Documents/1_university/grade_4/R/PES/allele_sample6_H1H2.csv")
fp_LR.mutation <- tclVar("C:/Users/Hide/Documents/1_university/grade_4/R/PES/GF.mutation_rates.csv")
Calculate(filepath_LR.frequency1,filepath_LR.DNAtype1,filepath_LR.mutation1,1,linkage.loci,RR)

  window <- tktoplevel()
  a <- tcltextvariable(window,value="awvawev")
  tkpack(a)
b <- tclVar("12345")
c <- "12356"
a <- tkgetSaveFile(parent = window, filetypes = "{{TXT Files} {.txt}}")
file.create(tclvalue(a))
write(paste(b,"\n",sep="",collapse=""),file=tclvalue(a))

paste(c("######DNA calculation result\n\n",
                      "######List of included systems\n",

                      "######List of involved persons\n"

                      ),sep=""
                    )
a <- c(1,2,3)
which(a == 0)
a[2] <- NA
length(a == 0)
b <- c("b","v","c")
paste(b,"\n",sep="",collapse="")

paste(a,collapse="",sep="")
    #tkmessageBox(message = LRprod)
  #label_result1 <- tklabel(frame1_result,textvariable=log_result1,font="our_font",width="30",anchor="w",justify="left",relief="groove",borderwidth="2")

a <- list()
a[[1]] <- list(2,3)
a[[2]] <- list(1,2)
sapply(a,function(x){return(x[[1]] + x[[2]])})
    pre_SaveFileName <- tkgetSaveFile(parent = window, filetypes = "{{PDF ファイル} {.pdf}}")
pre_SaveFileName <- paste(pre_SaveFileName,".pdf")
a <- "hogehoge"
pdf(file=pre_SaveFileName)
        SaveFileName <- tclvalue(pre_SaveFileName)

    if (is.element("RR", names(simulation_result))){
      log_Hp[[1]] <- c(log_Hp[[1]], paste(log_linkage.loci[,1], " - ", log_linkage.loci[,2]))
      log_Hp[[2]] <- c(log_Hp[[2]], signif(apply(LR.linkage.Hp,2,mean), digit=3))
      log_Hp[[3]] <- c(log_Hp[[3]], signif(apply(LR.linkage.Hp,2,median), digit=3))
      log_Hp[[4]] <- c(log_Hp[[4]], signif(apply(LR.linkage.Hp,2,max), digit=3))
      log_Hp[[5]] <- c(log_Hp[[5]], signif(apply(LR.linkage.Hp,2,min), digit=3))
      log_Hp[[6]] <- c(log_Hp[[6]], signif(apply(LR.linkage.Hp,2,sd), digit=3))

      log_Hd[[1]] <- c(log_Hd[[1]], paste(log_linkage.loci[,1], " - ", log_linkage.loci[,2]))
      log_Hd[[2]] <- c(log_Hd[[2]], signif(apply(LR.linkage.Hd,2,mean), digit=3))
      log_Hd[[3]] <- c(log_Hd[[3]], signif(apply(LR.linkage.Hd,2,median), digit=3))
      log_Hd[[4]] <- c(log_Hd[[4]], signif(apply(LR.linkage.Hd,2,max), digit=3))
      log_Hd[[5]] <- c(log_Hd[[5]], signif(apply(LR.linkage.Hd,2,min), digit=3))
      log_Hd[[6]] <- c(log_Hd[[6]], signif(apply(LR.linkage.Hd,2,sd), digit=3))
    }

      LR.linkage.Hp <- simulation_result[[LR.Hp]][[2]]
      LR.linkage.Hd <- simulation_result[[LR.Hd]][[2]]

 graph.draw1 <- function(){
#par(fig=c(0,0,0,0))
#par(mfrow=c(1,2)) 
#par(mai = c(0,0,0,0))
#par(omi=c(0,0,0,0))
#par(plt=c(0,1,0,1))
                                                          plot(1) 
                                                         }
                                          window.graph <- tktoplevel()
                                          frame <- tkframe(window.graph,bg="#00ff00",borderwidth="3",relief="groove")
                                          simulation_graph <- tkrplot(frame,graph.draw1, hscale=2, vscale=2)
                                          #frame1 <- tkframe(window.graph,bg="#00ffff",borderwidth="3",relief="groove")
                                          #label <- tklabel(frame1,text="text")
                                          tkgrid(frame)
                                         #tkpack(frame1)
                                          #tkpack(label)
                                          tkgrid(simulation_graph,padx="20")



    report_common <- paste(c("======== Input files ========\n",
                        "Frequency file :,",report_filepath_LR.frequency,"\n",
                        "Allele file :,",report_filepath_LR.DNAtype,"\n",
                        "Mutation file :,",report_filepath_LR.mutation,"\n",
                        "\n======== Locus set ========\n",
                        "Kit :,",report_locus_set_name,"\n",
                        "Locus set :,",report_allelenames,"\n",
                        "\n======== Hypothesis ========\n",
                        "Name,Fa,Ma\n",
                        report_pedigree
                      ),sep="",collapse=""
                    )
    report_LR.specific <- paste(c("\n\n======== Likelihood ========\n",
                              "総合尤度比 : ",report_LRprod,"\n\n",
                              report_result,"\n\n", 
                              report_linkage,"\n",
                              "\n======== Version ========\n",
                              "LRめいかー ver 0.0"
                            ),sep="",collapse=""
                          )
    report_simulation.specific <- paste(c("\n\n======== Likelihood ========\n",
                              "Hypothesis 1 ",report_Hp.title,"\n",
                              paste(report_Hp, collapse="\n"),
                              "\n\nHypothesis 2 ",report_Hd.title,"\n",
                              paste(report_Hd, collapse="\n"),
                              "\n======== Parameters ========\n",
                              "Mutation range : ",result[["range"]],"\n",
                              "\n\n======== Version ========\n",
                              "LRめいかー ver 0.0"
                            ),sep="",collapse=""
                          )

img <- tclVar()
a <- tkimage.create("photo",img,file="C:/Users/Hide/Documents/1_university/grade_4/R/GUI/chinos.png")
window <- tktoplevel()
f <- tklabel(window,image=a)
tkpack(f)

fc.H1 <- familycheck(famid=rep(1, length(df.H1[,"id"])), id = df.H1[,"id"], father.id = df.H1[,"pat"], mother.id = df.H1[,"mat"])
