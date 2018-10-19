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

function simPheno(popGeno,h2_1,h2_2,meanMaf,dist,parms,q1QTLs,q2QTLs,q12QTLs)
    @printf("read %.0f individuals and %.0f genotypes \n", size(popGeno,1),size(popGeno,2)-1)
    totQTLs = q1QTLs + q2QTLs + q12QTLs
    
    selectedLoci = []
    p = mean(convert(Array,popGeno[:,2:end]),1)/2.0
    while length(selectedLoci) < totQTLs
        oneLoci = sample(2:size(popGeno,2), 1, replace=false) #column of loci
        if in(oneLoci,selectedLoci) != true
            uniLoci = rand(Uniform(0,meanMaf))
            if(meanMaf-uniLoci)< p[oneLoci-1][] <= (meanMaf+uniLoci)  #this -1 is because p has 1 less length bec. popGeno has ID
                @printf("loci %.0f lower %.3f maf %.3f upper %.3f \n", oneLoci[]-1,meanMaf-uniLoci,p[oneLoci-1][],meanMaf+uniLoci)
                push!(selectedLoci,oneLoci)
            end
        end
    end
    @printf("mean MAF of selected loci: %.2f \n", mean(p[vcat(selectedLoci...).-1]))
    
    QTLs = vcat(selectedLoci...)   #columns of QTL since 1st column is ID
        
#    alpha = rand(Normal(0.0,1.0),totQTLs)
#    alpha = rand(Gamma(0.4,1.66),totQTLs)
    alpha = eval(parse("rand($dist$parms,$totQTLs)"))
    
    ###if you want zero mean BV
    Xc     = convert(Array{Float64},popGeno)
    Xc[:,2:end]   .-= ones(Float64,size(popGeno,1))*2p
    Qc     = Xc[:,QTLs]
    #### otherwise, Qc = convert(Array,popGeno[:,QTLs])
    
    u1 = Qc*(vcat([ones(q1QTLs), ones(q12QTLs), zeros(q2QTLs)]...).*alpha)
    vare1 = cov(u1)*(1-h2_1)/h2_1
    
    u2 = Qc*(vcat([zeros(q1QTLs), ones(q12QTLs), ones(q2QTLs)]...).*alpha)
    vare2 = cov(u2)*(1-h2_2)/h2_2
    
    e = rand(MvNormal([0.0; 0.0],[vare1 0;0 vare2]),size(popGeno,1))'

    y1 = 100 + u1 .+ e[:,1]
    y2 = 200 + u2 .+ e[:,2]
 
    G = cov([u1 u2])
    R = cov(e)
    h2sim  = Diagonal(G./(G+R))
    h2sim  = h2sim[find(h2sim)]
    
    println("genetic cov: $G")
    println("residual cov: $R")
    println("heritabilities: $h2sim")
    
    infoSimQTL = DataFrame(Any,size(popGeno,2)-1,6)
    snpID = [string(names(popGeno)[i]) for i in 2:size(popGeno,2)] #first one is ID
    infoSimQTL[:,1] .= snpID
    infoSimQTL[:,2] .= collect(1:length(snpID))
    infoSimQTL[:,3:end] = 0
    infoSimQTL[QTLs[1:(q1QTLs+q12QTLs)]-1,3] = 1
    infoSimQTL[QTLs[((q1QTLs+q12QTLs)+1):end]-1,4]   = 1
    infoSimQTL[QTLs[1:end]-1,5] = alpha #may need transpose
    infoSimQTL[:,6] .= vcat(p...)
    println(infoSimQTL)
    writecsv("infoSimQTL",convert(Array,infoSimQTL))
    phenoData = DataFrame(ID = Int64[], pheno1 = Float64[], pheno2 = Float64[], u1 = Float64[], u2 = Float64[], e1 = Float64[], e2 = Float64[])
    [push!(phenoData, [popGeno[row,:ID] y1[row] y2[row] u1[row] u2[row] e[row,:1] e[row,:2]]) for row in 1:length(y1)]
    writecsv("infoPheno",convert(Array,phenoData))
    @printf("returning phenotypes of %.0f individuals \n", size(phenoData,1))
    return snpID[find(infoSimQTL[:,5])], G, R, phenoData  #QTLs[1:end]-1 instead of returning QTL id
end

