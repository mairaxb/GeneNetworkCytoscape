rm(list=ls()) 
#--

data <-read.table("P_EXPR_antiobios.txt", header = T)
View(data)
#on retire la premiere colonne et on renomme les autres colonnes
data = data[,-1]
colnames(data)[1]="id"
colnames(data)[2]="H202_1"
colnames(data)[3]="H202_2"
colnames(data)[4]="H202_3"
colnames(data)[5]="amp_1"
colnames(data)[6]="amp_2"
colnames(data)[7]="amp_3"
colnames(data)[8]="gent_1"
colnames(data)[9]="gent_2"
colnames(data)[10]="gent_3"
colnames(data)[11]="kan_1"
colnames(data)[12]="kan_2"
colnames(data)[13]="kan_3"
colnames(data)[14]="norf_1"
colnames(data)[15]="norf_2"
colnames(data)[16]="norf_3"
colnames(data)[17]="none_1"
colnames(data)[18]="none_2"
colnames(data)[19]="none_3"
View(data)

#Lecture de la liste des gènes d’intérêt et chargement dans un vecteur
focus<-read.table("GOI.txt")
colnames(focus)[1]="id"
View(focus)
list_focus <-as.vector(focus[,"id"])
list_focus

#Filtrage des données d’expression correspondant aux locus d'intérêt
filter_list <- data[,"id"] %in% list_focus #si il y a un match = True sinon False
data_focus<-data[filter_list,] #on ne va garder que les match donc les True
data_focus
#on ne garde que les locus d'intérets (locus connus pour être resistants aux ATB)

#--
#Filtrage des donnees pour n'avoir que les genes impliques dans utilisation du fer
#on veut etudier relation entre utilisation du fer et resistance aux ATB

fur1 <-read.table("GOI_FUR_1.txt", header = T)
fur2 <-read.table("GOI_FUR_2.txt", header = T)
View(fur1)
View(fur2)

fur1 = fur1[,-1]
fur1

fur2 = fur2[,-1]
fur2

#on ne garde que les genes avec sous/sur expression de 2
fur1=fur1[abs(fur1["log2.fold_change._Dfur.wt."])>2,]
fur2=fur2[abs(fur2["log2.fold_change._.Dfur.wt."])>2,]
View(fur1)
View(fur2)

#
colnames(fur1)[1]="id"
colnames(fur2)[1]="id" 

vect_fer <-as.vector(fur1[,"id"])
vect_fer

vect_atb <-as.vector(fur2[,"id"])
vect_atb

filter_1 <- data[,"id"] %in% vect_fer #si il y a un match = True sinon False
data_fer<-data[filter_1,] #on ne va garder que les match donc les True
data_fer
View(data_fer)

filter_2 <- data[,"id"] %in% vect_atb #si il y a un match = True sinon False
data_atb<-data[filter_2,] #on ne va garder que les match donc les True
data_atb
View(data_atb)

#--
#tout mettre dans un seul fichier pour faire une correlation 
#elimination des doublons
new_l=rbind(data_focus, data_fer, data_atb)
doublons <- which(duplicated(new_l$id)) 
new_l<-new_l[-doublons,] 
View(new_l)

new_l2 = new_l
rownames(new_l2)=new_l[,1]
new_l2=new_l2[,-1]
new_l2=t(new_l2)
View(new_l2)

result_correlation=cor(new_l2)
View(result_correlation)

#construction de la matrice de correlation
library("reshape2")
list_correlation=melt(result_correlation)
list_correlation
View(list_correlation)
list_correlation=list_correlation[abs(list_correlation["value"])>0.93,]
#Il est nécessaire de supprimer les paires redondantes (on est partie d'une matrice symétrique)
#Création d'une nouvelle colonne Alphabétique pour stipuler si le nome du gene A est ordonné selon l'aphabet par rapport au nom du gene B
list_correlation["Alphabétique"]<-as.character(list_correlation[,"Var1"])<as.character(list_correlation[,"Var2"])
View(list_correlation)

#Supression des lignes où l'ordre de classement entre les deux noms de gènes ne suit par l'ordre alphabétique (permet de supprimer les doublons)
list_correlation=list_correlation[list_correlation[,4]==TRUE,]
list_correlation
View(list_correlation)

#Supression de la colonne transitoire alphabétique - n'est plus nécessaire-
list_correlation=list_correlation[,-4]
#Supression de la colonne Value
list_correlation=list_correlation[,-3]
#Ajout de la colonne Interaction pour donner le type, ici ce sont les paires dont  l'expression est corrélée.
list_correlation["type"]="correlation_pair"
list_correlation
LISTE = list_correlation 

list_interaction <-read.table("sRNA_interaction.txt", header = T)
View(list_interaction)
vector_coli <-as.vector(new_l[,"id"])
vector_coli

filter_coli <- list_interaction[,"ID"] %in% vector_coli #si il y a un match = True sinon False
data_RNA <- list_interaction[filter_coli,] #on ne va garder que les match donc les True
data_RNA

colnames(data_RNA)[1]="Var1"
colnames(data_RNA)[2]="Var2"
colnames(data_RNA)[3]="type"
data_RNA

#Fusion des deux listes pour créer un nouveau data frame et le sauvgarder dans un fichier
FINAL=rbind(LISTE,data_RNA)
FINAL

write.table (FINAL,"FINAL.txt", sep="\t",row.names=FALSE,quote=F)

#--
FT <-read.table("P_TF_LABEL_COLI.txt", header = T)
colnames(FT)[1]="genes"
FT

list_coli <-read.table("sRNA_list_coli.txt", header = T)
list_coli 
colnames(list_coli)[1]="genes"

#Fusion des deux listes pour créer un nouveau data frame et le sauvgarder dans un fichier
REG=rbind(FT,list_coli)
REG
write.table (REG,"REG.txt", sep="\t",row.names=FALSE,quote=F)
