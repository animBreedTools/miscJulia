using Distributions
#using RCall

function samplePop(genotypes,whichGen,snpInfo,chr,nRef,nTest)
    @printf("read %.0f individuals and %.0f genotypes \n", size(genotypes,1),size(genotypes,2)-1)
    tempMapData = readtable(snpInfo,header=false,separator=' ')
    tempMapData = tempMapData[tempMapData[:x4].<=chr,:]
    nTot = nRef+nTest
    if length(whichGen)==1
        genLims = ((2200*whichGen[1])-2200+1):(2200*whichGen[1])
        println("these individuals were used to select populations ",genLims)
        allInd  = sample(genLims, nTot, replace=false)
        refInd  = allInd[1:nRef]
        testInd = allInd[(nRef+1):end]
    end
    if length(whichGen)>1
    refGener   = 1:(whichGen[1,:][end]*2200)
    testGener  = ((whichGen[1,:][end]*2200)+1):(whichGen[2,:][end]*2200)
    refInd     = sample(refGener, nRef, replace=false)
    testInd    = sample(testGener, nTest, replace=false)
    allInd     = [refInd ;testInd]
    end
    popGeno = genotypes[allInd,1:(size(tempMapData,1)+1)]
    @printf("returning %.0f individuals and %.0f genotypes on %.0f chromosomes \n", size(popGeno,1),size(popGeno,2)-1,chr)
    return nTot, refInd, testInd, popGeno
end 

function simPheno(popGeno,h2_1,h2_2,meanMaf,dist,parms,q1QTLs,q2QTLs,q12QTLs,snpInfo,nChr)
    @printf("read %.0f individuals and %.0f genotypes \n", size(popGeno,1),size(popGeno,2)-1)
    totQTLs = q1QTLs + q2QTLs + q12QTLs
    
    mapData = readtable(myMap,header=false,separator=' ')
    println(by(mapData, :x4, nrow))
    chrData = cumsum(by(mapData, :x4, nrow)[:x1])[nChr]
    
    usedMapData = mapData[1:chrData,:]
    
    finalGroups = []
    count = 0

for thisChr in 1:nChr
    thisRange1 = collect(0:1000000:(maximum(usedMapData[usedMapData[:x4].==thisChr,:x5])))
    thisRange2 = collect(1000000:1000000:(maximum(usedMapData[usedMapData[:x4].==thisChr,:x5])+1000000))
    finalRange = [thisRange1 thisRange2]

    theseSNPs = usedMapData[usedMapData[:x4].==thisChr,:]
for groupChr in 1:size(finalRange,1)
    groupSize = length(theseSNPs[finalRange[groupChr,1].<theseSNPs[:x5].<=finalRange[groupChr,2],:x5])
    if groupSize !=0
            count = count+1
            finalGroups = push!(finalGroups,vcat(fill(count,groupSize))...)
    end
end    
end
    
    usedMapData.groups = finalGroups
    usedMapData.nSNP   = collect(1:size(usedMapData,1))
    
    selectedLoci = []
    selectedBins = []

    p = mean(convert(Array,popGeno[1:2200,2:end]),1)./2  ##### #only based IND 
    q = 1-p
    minPQ = copy(p)
    for i in 1:length(p)
        minPQ[i] = min(p[i],q[i])
    end
    
    while length(selectedLoci) < totQTLs
        thisBin = sample(usedMapData[:groups],1)
        println(usedMapData[usedMapData[:groups].==thisBin[],:])
        if in(thisBin,selectedBins) != true
            oneLoci = sample(usedMapData[usedMapData[:groups].==thisBin[],:nSNP], 1)
            if in(oneLoci,selectedLoci) != true
                uniLoci = rand(Uniform(0,meanMaf))   #0.01<maf<0.30 I put 0.05 instead of "meanMaf"
                if(meanMaf-uniLoci)< minPQ[oneLoci][] <= (meanMaf+uniLoci)
                    @printf("bin %.0f loci %.0f lower %.3f maf %.3f upper %.3f \n", thisBin[], oneLoci[],meanMaf-uniLoci,minPQ[oneLoci][],meanMaf+uniLoci)
                    push!(selectedLoci,oneLoci)
                    push!(selectedBins,thisBin[])
                end
            end
    end
end
    @printf("mean MAF of selected loci: %.2f \n", mean(minPQ[vcat(selectedLoci...)]))
QTLs = shuffle(vcat(selectedLoci...))
 
##(a) alpha = eval(parse("rand($dist$parms,$totQTLs)"))
    
   # alpha = rcopy(R"""
   # library(lcmix)
   # rho <- matrix(c(1,0.9,0.9,1),nrow=2)
   # alpha <- rmvgamma($totQTL,0.4,1.66,rho)
   # """)
    
    alpha = convert(Array,readtable("rSimAlpha",header=false,separator=' '))

    nPosCor = Int(ceil(q12QTLs*0.78))
    nNegCor = q12QTLs - nPosCor
    dNeg1   = sample([-1 1],nNegCor)
    dNeg2   = -1*dNeg1
    dPos    = sample([-1 1],nPosCor)
    dNegCor = [dNeg1 dNeg2]
    dPosCor = [dPos dPos]
    dNegPos = [dNegCor;dPosCor]
    dNegPos = dNegPos[shuffle(1:end), :]

    alpha = alpha.*[sample([-1 1],q1QTLs) zeros(q1QTLs);dNegPos;zeros(q2QTLs) sample([-1 1],q2QTLs)] 
       
    Xc     = convert(Array{Float64},popGeno)    
    Xc[:,2:end]   .-= ones(Float64,size(popGeno,1))*2p
    Qc     = Xc[:,QTLs]
    
    u = Qc*alpha
    
    u1 = u[:,1]
    u2 = u[:,2]
    
    #----------------
    Gnow  = cov([u1[1:2200] u2[1:2200]])
    println("alpha $(cov(alpha))")
    println("Gnow $Gnow")
    c     = diag(100.0 ./Gnow)
    alpha = sqrt.(c').*alpha
    u = Qc*alpha
    u1 = u[:,1]
    u2 = u[:,2]
    println("scaled alpha $(cov(alpha))")
    println("scaled G $(cov([u1[1:2200] u2[1:2200]]))")
    #----------------
    
    vare1 = var(u1[1:2200])*(1-h2_1)/h2_1    #only based IND
    vare2 = var(u2[1:2200])*(1-h2_2)/h2_2
    
 ##(a)   u1 = Qc*(vcat([ones(q1QTLs), ones(q12QTLs), zeros(q2QTLs)]...).*alpha)
 ##(a)   vare1 = cov(u1)*(1-h2_1)/h2_1
 ##(a)   
 ##(a)   u2 = Qc*(vcat([zeros(q1QTLs), ones(q12QTLs), ones(q2QTLs)]...).*alpha)
 ##(a)   vare2 = cov(u2)*(1-h2_2)/h2_2
    
    e = rand(MvNormal([0.0; 0.0],[vare1 0;0 vare2]),size(popGeno,1))'
        
    y1 = u1 .+ e[:,1]                       #no add Mean
    y2 = u2 .+ e[:,2]
 
    G = cov([u1[1:2200] u2[1:2200]])      #only based IND
    R = cov(e[1:2200,:])                  #only based IND
    h2sim  = Diagonal(G./(G+R))
    h2sim  = h2sim[find(h2sim)]
    
    println("genetic cov: $G")
    println("residual cov: $R")
    println("heritabilities: $h2sim")
    
    infoSimQTL = DataFrame(Any,size(popGeno,2)-1,7)
    snpID = [string(names(popGeno)[i]) for i in 2:size(popGeno,2)] #first one is ID
    infoSimQTL[:,1] = snpID
    infoSimQTL[:,2] = collect(1:length(snpID))
    infoSimQTL[:,3:end] = 0
    
    infoSimQTL[QTLs[1:(q1QTLs+q12QTLs)],3] = 1000         #wrong because the order of QTLs change up
    infoSimQTL[QTLs[(q1QTLs+1):end],4] = 1000   #wrong because the order of QTLs change up
    infoSimQTL[QTLs[1:end],5:6] = DataFrame(alpha)

    infoSimQTL[:,7] = vcat(p...)         #this was 6
    println(infoSimQTL)
    writecsv("infoSimQTL",convert(Array,infoSimQTL))
    phenoData = DataFrame(ID = Int64[], pheno1 = Float64[], pheno2 = Float64[], u1 = Float64[], u2 = Float64[], e1 = Float64[], e2 = Float64[])
    [push!(phenoData, [popGeno[row,:ID] y1[row] y2[row] u1[row] u2[row] e[row,:1] e[row,:2]]) for row in 1:length(y1)]
    writecsv("infoPheno",convert(Array,phenoData))
    @printf("returning phenotypes of %.0f individuals \n", size(phenoData,1))
    return snpID[unique(vcat(find(infoSimQTL[:,5]),find(infoSimQTL[:,6])))], G, R, phenoData  #QTLs[1:end]-1 instead of returning QTL id
end

