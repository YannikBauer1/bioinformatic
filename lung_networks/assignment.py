import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import pylab
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


#----------------------------------------------------------#
#---------------------- Task 1 ----------------------------#
#----------------------------------------------------------#

# os prints e graficos sao todos em commentarios para nao ser demaisado em uma vez

class myClass:
    def __init__(self, path):
        self.path = path
        self.df = self.createDf(path)
        self.datatype = path.split("_")[1]
        self.corr = ["pearson","spearman"]
        self.corrDic = {"pearson":0,"spearman":1}
        self.corrMatriz = []
        self.X = [0.1,0.2,0.3,0.4,0.5]
        self.graphs = []

    def createDf(self,path):
        df = pd.read_csv(path, delimiter="\t")
        new_col = list(df.columns.values[1:])
        new_col = [i[:12] for i in new_col] + ["out"]
        df = df.set_axis(new_col, axis=1, inplace=False)
        df = df.drop(["out"], axis=1)
        return df

    def createCorrMatriz(self):
        l=[]
        for i in self.corr:
            l.append(self.df.corr(method=i))
        self.corrMatriz=l

    def printCorrMatriz(self,corr=""):
        if self.corrMatriz==[]:
            self.createCorrMatriz()
        correlations = []
        if type(corr) == list:
            correlations = corr
        elif corr=="":
            correlations = self.corr
        else:
            correlations.append(corr)
        for i in correlations:
            if i == "pearson":
                print("Correlation Matriz | datatype: ",self.datatype," | correlation=",i," | :")
                print(self.corrMatriz[0])
            if i == "spearman":
                print("Correlation Matriz | datatype: ",self.datatype," | correlation=",i," | :")
                print(self.corrMatriz[1])

    def createThresholdMatriz(self,new_df,x):
        df = new_df.copy()
        columns = list(df.columns.values)
        for col in columns:
            df.loc[df[col] >= x, col] = 1
            df.loc[df[col] <= -x, col] = 1
            df.loc[df[col] != 1, col] = 0
            df.at[col, col] = 0
        return df

    def printThresholdMatriz(self,x,corr=""):
        if self.corrMatriz==[]:
            self.createCorrMatriz()
        correlations = []
        if type(corr) == list:
            correlations = corr
        elif corr=="":
            correlations = self.corr
        else:
            correlations.append(corr)
        for i in correlations:
            if i == "pearson":
                print("Correlation Matriz with threshold ",x," | datatype: ",self.datatype," | correlation=",i," | :")
                print(self.createThresholdMatriz(self.corrMatriz[0],x))
            if i == "spearman":
                print("Correlation Matriz with threshold ",x," | datatype: ",self.datatype," | correlation=",i," | :")
                print(self.createThresholdMatriz(self.corrMatriz[1],x))

    def createGraph(self,x,corr=""):
        if self.corrMatriz==[]:
            self.createCorrMatriz()
        correlations = []
        if type(corr) == list:
            correlations = corr
        elif corr=="":
            correlations = self.corr
        else:
            correlations.append(corr)

        graphs = []
        for i in correlations:
            index = self.corrDic[i]
            m = self.createThresholdMatriz(self.corrMatriz[index],x)
            graph = nx.from_numpy_matrix(m.values)
            columns = list(self.df.columns.values)
            mapping = {i: columns[i] for i in range(len(columns))}
            graph = nx.relabel_nodes(graph, mapping)
            graphs.append(graph)
        return graphs

    def printGraph(self,x,corr=""):
        graphs = self.createGraph(x, corr)
        correlations = []
        if type(corr) == list:
            correlations = corr
        elif corr=="":
            correlations = self.corr
        else:
            correlations.append(corr)
        for i in range(len(correlations)):
            nx.draw(graphs[i], with_labels=True)
            fig = pylab.gcf()
            fig.canvas.manager.set_window_title(self.datatype+" | "+correlations[i]+" | "+str(x))
            plt.show()

    def createGraphs(self):
        graphs=[]
        for x in self.X:
            g = self.createGraph(x)
            graphs.append(g)
        self.graphs=graphs

    def printStatisticsTable(self):
        if self.graphs == []:
            self.createGraphs()
        for index in list(self.corrDic.values()):
            for g in range(len(self.graphs)):
                graph = self.graphs[g][index]
                x = self.X[g]
                n_nodes = graph.number_of_nodes()
                n_edges = graph.number_of_edges()
                avg_degree = sum([k[1] for k in graph.degree()]) / n_nodes
                avg_cc = nx.average_clustering(graph)
                print("{:<3} {:<8} {:<4} {:<4} {:>.2f} {:>.2f}".format(x, self.corr[index], n_nodes, n_edges, avg_degree,avg_cc))

    def networkAnalytics(self,x,corr=""):
        graphs = self.createGraph(x, corr)
        correlations = []
        if type(corr) == list:
            correlations = corr
        elif corr == "":
            correlations = self.corr
        else:
            correlations.append(corr)

        for i in range(len(correlations)):
            c=correlations[i]
            g=graphs[i]
            p_k = nx.degree_histogram(g)
            p_k = [round(k / sum(p_k), 3) for k in p_k]
            indexes = [k for k in range(len(p_k))]
            p_c = []
            for y in indexes:
                nodes_degree_k = [u for u in g.nodes() if g.degree(u) == y]
                clustering_coeff = nx.clustering(g, nodes_degree_k)
                if len(clustering_coeff) == 0:
                    p_c.append(0)
                else:
                    p_c.append(sum(clustering_coeff.values()) / len(clustering_coeff))
            fig = pylab.gcf()
            fig.canvas.manager.set_window_title(self.datatype + " | " + correlations[i] + " | " + str(x))
            plt.subplot(1,2,1)
            plt.title("P(k) for each degree")
            plt.scatter(indexes, p_k)
            plt.subplot(1,2,2)
            plt.title("C(k) for each degree")
            plt.scatter(indexes, p_c)
            plt.show()

gene = myClass("Lung/LUNG_Gene_Expression.txt")
mirna = myClass("Lung/LUNG_Mirna_Expression.txt")
methy = myClass("Lung/LUNG_Methy_Expression.txt")

### task 1.1
#gene.printCorrMatriz()
#mirna.printCorrMatriz()
#methy.printCorrMatriz()

### task 1.2
#gene.printGraph(0.3)
#mirna.printGraph(0.3)
#methy.printGraph(0.3)

### task 1.3
#gene.printCorrMatriz(corr="pearson")
#gene.printThresholdMatriz(0.3,corr="pearson")
#gene.printCorrMatriz(corr="spearman")
#gene.printThresholdMatriz(0.3,corr="spearman")
#mirna.printCorrMatriz(corr="pearson")
#mirna.printThresholdMatriz(0.3,corr="pearson")
#mirna.printCorrMatriz(corr="spearman")
#mirna.printThresholdMatriz(0.3,corr="spearman")
#methy.printCorrMatriz(corr="pearson")
#methy.printThresholdMatriz(0.3,corr="pearson")
#methy.printCorrMatriz(corr="spearman")
#methy.printThresholdMatriz(0.3,corr="spearman")

### task 1.4
#gene.printStatisticsTable()
#mirna.printStatisticsTable()
#methy.printStatisticsTable()

### task 1.5
#gene.networkAnalytics(0.2)
#gene.printGraph(0.2)
#print("Olhando para os graficos de 'gene' o network tipo mais semelhante parece ser Scale-free Network")

#mirna.networkAnalytics(0.2)
#mirna.printGraph(0.2)
#print("Olhando para os graficos de 'mirna' o network tipo mais semelhante parece ser Scale-free Network")

#methy.networkAnalytics(0.2)
#methy.printGraph(0.2)
#print("Olhando para os graficos de 'methy' é mais dificil, no grafico que mostra o C(k) para cada degree\n"
#      "parece ser semelhante ao  hierarchical, olhando para o grafico que mostra o P(k) para cada degree\n"
#      "e para o grafo em si podia também ser scale-free ou random")



#----------------------------------------------------------#
#---------------------- Task 2 ----------------------------#
#----------------------------------------------------------#

# a chamada da funcao task2 está em commentario, para poder chamar quando é preciso
# e nao mostrar logo tudo

def elbowMethod(all_data):
    Sum_of_squared_distances = []
    K = range(1, 10)
    for num_clusters in K:
        kmeans = KMeans(n_clusters=num_clusters)
        kmeans.fit(all_data)
        Sum_of_squared_distances.append(kmeans.inertia_)
    plt.plot(K, Sum_of_squared_distances, "bx-")
    plt.xlabel("Values of K")
    plt.ylabel("Sum of squared distances / Inertia")
    plt.title("Elbow Method For Optimal k")
    plt.show()

def silhouetteAnalysis(all_data):
    K = range(2, 10)
    silhouette_avg = []
    for num_clusters in K:
        kmeans = KMeans(n_clusters=num_clusters)
        kmeans.fit(all_data)
        cluster_labels = kmeans.labels_
        silhouette_avg.append(silhouette_score(all_data, cluster_labels))
    plt.plot(K, silhouette_avg, "bx-")
    plt.xlabel("Values of K")
    plt.ylabel("Silhouette score")
    plt.title("Silhouette analysis For Optimal k")
    plt.show()

def task2():
    df1 = gene.df
    df2 = mirna.df
    df3 = methy.df

    all_data = pd.concat([df1,df2,df3])
    all_data = all_data.transpose()

    #elbowMethod(all_data)
    # k 3 é o melhor

    #silhouetteAnalysis(all_data)
    # k 2 é o melhor e depois 3 mas com silhouette score baixo

    model = KMeans(n_clusters=3)
    model.fit(all_data)
    values=list(model.labels_)

    all_data_corr = all_data.transpose().corr()
    all_data_corr01 = all_data_corr.copy()
    columns = list(all_data_corr01.columns.values)
    for col in columns:
        all_data_corr01.loc[all_data_corr01[col] >= 0.2, col] = 1
        all_data_corr01.loc[all_data_corr01[col] <= -0.2, col] = 1
        all_data_corr01.loc[all_data_corr01[col] != 1, col] = 0
        all_data_corr01.at[col, col] = 0
    graph = nx.from_numpy_matrix(all_data_corr01.values)
    columns = list(all_data_corr01.columns.values)
    mapping = {i: columns[i] for i in range(len(columns))}
    graph = nx.relabel_nodes(graph, mapping)
    nx.draw(graph,cmap=plt.get_cmap('viridis'), node_color=values, with_labels=True)
    plt.show()

    df4 = pd.read_csv("Lung/LUNG_Survival.txt", delimiter="\t")
    df4 = df4.transpose()
    coln = list(df4.loc[["PatientID"]].values)[0]
    coln = [i[:12] for i in coln]
    df4 = df4.set_axis(coln, axis=1, inplace=False)
    df4 = df4.drop("PatientID")
    df5 = pd.read_csv("Lung/lusc.clinical.txt", delimiter="\t")
    df5 = df5.transpose()
    coln = list(df5.loc[["sample"]].values)[0]
    coln = [i[:12].replace(".", "-") for i in coln]
    df5 = df5.set_axis(coln, axis=1, inplace=False)
    df5 = df5.drop("sample")
    aux_data = pd.concat([df4, df5])

    cluster1=[columns[i] for i in range(len(columns)) if values[i]==0]
    cluster2=[columns[i] for i in range(len(columns)) if values[i]==1]
    cluster3=[columns[i] for i in range(len(columns)) if values[i]==2]

    aux_data1=aux_data[cluster1].transpose()
    aux_data2=aux_data[cluster2].transpose()
    aux_data3=aux_data[cluster3].transpose()

    a1 = aux_data1.age_at_diagnosis.astype(float).describe()
    a2 = aux_data2.age_at_diagnosis.astype(float).describe()
    a3 = aux_data3.age_at_diagnosis.astype(float).describe()
    print(pd.concat([a1,a2,a3],axis=1).set_axis(["cluster1","cluster2","cluster3"], axis=1, inplace=False))
    # nao há muito diferenca entre os valores

    a1 = aux_data1.Survival.astype(float).describe()
    a2 = aux_data2.Survival.astype(float).describe()
    a3 = aux_data3.Survival.astype(float).describe()
    print(pd.concat([a1,a2,a3],axis=1).set_axis(["cluster1","cluster2","cluster3"], axis=1, inplace=False))
    # há uma grande diferenca entre o mean, median q1, q3 para os clusters

    a1 = aux_data1.pathologic_stage.value_counts()
    a1index= list(a1.index)
    a1 = [i/sum(a1.values) for i in a1.values]
    a1Dict = {a1index[i]: a1[i] for i in range(len(a1index))}
    a2 = aux_data2.pathologic_stage.value_counts()
    a2index = list(a2.index)
    a2 = [i / sum(a2.values) for i in a2.values]
    a2Dict = {a2index[i]: a2[i] for i in range(len(a2index))}
    a3 = aux_data3.pathologic_stage.value_counts()
    a3index = list(a3.index)
    a3 = [i / sum(a3.values) for i in a3.values]
    a3Dict = {a3index[i]: a3[i] for i in range(len(a3index))}
    a0 = a1index+a2index+a3index
    a0 = list(dict.fromkeys(a0))
    dic = {i: [] for i in a0}
    for j in [a1Dict,a2Dict,a3Dict]:
        for k in a0:
            if k in j:
                dic[k]= dic[k]+[j[k]]
            else:
                dic[k]= dic[k]+[0]
    print(pd.DataFrame.from_dict(dic).transpose().set_axis(["cluster1","cluster2","cluster3"], axis=1, inplace=False))
#task2()



