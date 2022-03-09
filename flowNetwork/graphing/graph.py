def bokeh_interactive_matrix_plot(mat,nonzero_only=True,
                                    plotWidth=640,plotHeight=640,
                                    colorMap=sns.color_palette("coolwarm", n_colors=256).as_hex(),
                                    vmax=None,vmin=None,
                                    xlabel=None,ylabel=None,title=None
                                 ):
    if nonzero_only:
        matInds=np.nonzero(mat)
        matVals=np.array(mat[matInds]).flatten()
    else:
        matInds=(np.arange(mat.shape[0]),np.arange(mat.shape[1]))
        matVals=mat[matInds]
    
    tempFrame=pd.DataFrame({
        'X':matInds[0],
        'Y':matInds[1],
        'Val':matVals
    })
    
    vMax=np.max(matVals) if vmax is None else vmax
    vMin=np.min(matVals) if vmin is None else vmin
        
    valueColName='Val'
    
    mapper=LinearColorMapper(palette=colorMap,
                         low=vMin,
                         high=vMax)

    heatMapColumnSource=ColumnDataSource(tempFrame)

    toolTips=[(colName,"@{}".format(colName)) for colName in tempFrame.columns]

    p=bk.plotting.figure(plot_width=plotWidth,plot_height=plotHeight,
                tooltips=toolTips)

    if not (xlabel is None):
        p.xaxis.axis_label=xlabel
    if not (ylabel is None):
        p.yaxis.axis_label=ylabel
    if not (title is None):
        p.title.text=title
    
    p.background_fill_color="black"
    p.background_fill_alpha=.75
    
    color_bar=ColorBar(color_mapper=mapper,location='top_left')

    p.rect(x='X',y='Y',
           width=1,height=1,
           source=heatMapColumnSource,
           fill_color={'field':valueColName,'transform':mapper},line_color=None)
    p.add_layout(color_bar)
    show(p)
    
#A relatively simple bokeh based heatmap tool to quickly plot a heat map
#using a given value column and coordinate column pair from a data table
def bokeh_dataTable_heatMap(plotData,Xcol,Ycol,dataCol,
                            width=640,height=640,
                            rectheight=1,rectwidth=1,
                            colorMap=sns.color_palette("coolwarm", n_colors=256).as_hex(),
                            title=None,
                            xlabel=None,
                            ylabel=None,
                            axisFontSize="14pt",
                            vmin=None,vmax=None):
    
    if vmin is None:
        eMin=plotData[dataCol].min()
    else:
        eMin=vmin
    if vmax is None:
        eMax=plotData[dataCol].max()
    else:
        eMax=vmax
        
    plotData['ColorWeight']=plotData[dataCol] #.map(lambda x: cNorm(x)) #.map(np.log10)
    
    p=bk.plotting.figure(
        plot_width=width,plot_height=height)
    
    if title:
        p.add_layout(Title(text=title),'above')
    if xlabel:
        p.xaxis.axis_label=xlabel
    if ylabel:
        p.yaxis.axis_label=ylabel
        
    p.xaxis.major_label_text_font_size=axisFontSize
    p.yaxis.major_label_text_font_size=axisFontSize
    
    src=ColumnDataSource(plotData)
    mapper=linear_cmap(field_name='ColorWeight',
                       palette=colorMap,
                       low=eMin,
                       high=eMax,
                      )
    
    cbMapper=LinearColorMapper(palette=colorMap,
                               low=eMin,
                               high=eMax)
    color_bar = ColorBar(color_mapper=cbMapper, ticker= BasicTicker(),
                     location=(0,0))

    p.add_layout(color_bar, 'right')

    p.rect(source=src,width=rectwidth,height=rectheight, x=Xcol,y=Ycol,
             fill_color=mapper,color=mapper)

    
    
    hover=HoverTool()
    hover.tooltips=[
        (colName,"@{"+"{colName}".format(colName=colName)+"}")
        for colName in src.data.keys()
    ]
    print(hover.tooltips)
    p.add_tools(hover)

    output_notebook()
    show(p)
    
#functions for drawing networks using NGLview
import nglview as nv
import pytraj as pt

def bokeh_tidtyDataTable_barChart(dataTable,xColumns,valColumn,
                                  fillColumnInd=0,lineColumnInd=0,
                                  fillPalette=Spectral11,linePalette=Spectral11,
                                  columnFormatFuns=None,
                                  xAxisLabel=None,yAxisLabel=None,
                                  pWidth=500,pHeight=900):
    plotData=dataTable[np.concatenate([xColumns,[valColumn]])].copy()
    if columnFormatFuns:
        for iCol,xCol in enumerate(xColumns):
            plotData[xCol]=plotData[xCol].map(columnFormatFuns[iCol])
    else:
        for xCol in xColumns:
            plotData[xCol]=plotData[xCol].map(str)
            
    xVals=[tuple(x) for x in list(plotData[xColumns].set_index(xColumns[0]).to_records())]
    
    source=ColumnDataSource(data={
        'x':xVals,
        valColumn:plotData[valColumn]
    })
    
    ToolTips=[(colName,'@%s'%colName) for colName in \
              source.to_df().columns]
    
    p=figure(x_range=FactorRange(*xVals),plot_height=500,plot_width=900,tooltips=ToolTips)
    p.vbar(x='x',top=valColumn,width=.9,source=source,
           fill_color=factor_cmap('x',palette=fillPalette,
                                  factors=plotData[xColumns[fillColumnInd]].unique(),
                                  start=fillColumnInd,end=fillColumnInd+1),
           line_color=factor_cmap('x',palette=linePalette,
                                  factors=plotData[xColumns[lineColumnInd]].unique(),
                                  start=lineColumnInd,end=lineColumnInd+1))
    
    p.xaxis.major_label_orientation = 87.5*np.pi/180.
    p.xaxis.group_label_orientation = 85.*np.pi/180.
    p.xaxis.subgroup_label_orientation = 80.*np.pi/180.
    if xAxisLabel:
        p.xaxis.axis_label = xAxisLabel
    if yAxisLabel:
        p.yaxis.axis_label = yAxisLabel
    show(p)


def getCorrNetEdgeColors(valMat,maskInds=None,
                         cmap=matplotlib.cm.get_cmap('viridis',1024),
                         cNorm=None,
                         baseColor=[0.,0.,0.]):
    if maskInds is None:
        nzInds=np.nonzero(np.isfinite(valMat))
    else:
        nzInds=maskInds
    
    if cNorm is None:
        vMin=np.min(valMat[nzInds])
        vMax=np.max(valMat[nzInds])
        norm=matplotlib.colors.Normalize(vmin=vMin,vmax=vMax)
    else:
        norm=cNorm
    
    tempCarray=cmap(norm(valMat[nzInds]).data)
    edgeColors=np.zeros((valMat.shape[0],valMat.shape[1],3))
    for ii in np.arange(3):
        edgeColors[:,:,ii]=baseColor[ii]
        edgeColors[:,:,ii][nzInds]=tempCarray[:,ii]
    return(edgeColors)

def drawEdgeColorBar(valMat,maskInds=None,
                     cmap=matplotlib.cm.get_cmap('viridis',1024),
                     cNorm=None,
                     ax=None,
                     barSize=[16,4],
                     label='Edge Color Scale',
                     label_size=10,
                     nTicks=5,
                     orientation='horizontal',
                     verbose=False):
    if maskInds is None:
        nzInds=np.nonzero(np.isfinite(valMat))
    else:
        nzInds=maskInds
    tempCmap=cmap
    tempCarray=tempCmap(valMat[nzInds])
    if verbose:
        print('vMin:',vMin,'vMax:',vMax)
    
    if cNorm is None:
        vMin=np.min(valMat[nzInds])
        vMax=np.max(valMat[nzInds])
        norm=matplotlib.colors.Normalize(vmin=vMin,vmax=vMax)
    else:
        norm=cNorm
    
    if ax is None:
        if orientation=='horizontal':
            fig = plt.figure(figsize=barSize)
            ax1 = fig.add_axes([0.3, 0.80, 0.4, 0.15])
        else:
            fig = plt.figure(figsize=np.flip(barSize))
            ax1 = fig.add_axes([0.3, 0.10, 0.15, 0.4])
    else:
        ax1=ax
    cb1 = matplotlib.colorbar.ColorbarBase(ax1, cmap=tempCmap,
                                    norm=norm,
                                    orientation=orientation)
    cb1.set_label(label,size=label_size)
    cb1.ax.tick_params(labelsize=label_size)
    tick_locator = matplotlib.ticker.MaxNLocator(nbins=nTicks)
    cb1.locator = tick_locator
    cb1.update_ticks()
    plt.show()

def getCorrNetEdgeRadii(valMat,maskInds=None,
                        eMin=.0625,eMax=.75):
    if maskInds is None:
        nzInds=np.isfinite(valMat)
    else:
        nzInds=maskInds
        
    plotMat=valMat
    radiiMat=np.zeros(plotMat.shape)
    radiiMat[nzInds]=np.log(np.abs(np.array(plotMat)[nzInds]))
    rMin=np.min(radiiMat[nzInds])
    rMax=np.max(radiiMat[nzInds])
    rLinFun=lambda x: (x-rMin)/(rMax-rMin)
    rMap=lambda x:(eMax-eMin)*rLinFun(x)+eMin
    radiiMat[nzInds]=rMap(np.array(radiiMat[nzInds]))
    return(radiiMat)

def drawProtCorrMat(protStruc,corrMat,ngViewOb,
                    frame=0,colorsArray=None,radiiMat=None,
                    undirected=True):
    nzInds=np.nonzero(corrMat)
    edgeList=[]
    if not (radiiMat is None):
        radMat=radiiMat
    else:
        radMat=(np.abs(corrMat)>0)*.25
    if undirected:
        indArray=np.array([[ind[0],ind[1]] \
                            for ind in np.array([nzInds[0],nzInds[1]]).T \
                            if ind[0]<ind[1]])
    else:
        indArray=np.array([nzInds[0],nzInds[1]]).T
    for nzInd in indArray:
        if not (colorsArray is None):
            colorParm={'edgeColor':list(colorsArray[nzInd[0],nzInd[1],:])}
        else:
            colorParm={}
        edgeList.append(drawProtNetEdge(
            protStruc,nzInd[0]+1,nzInd[1]+1,
            ngViewOb,frame,radius=radMat[nzInd[0],nzInd[1]],
            **colorParm))
    return edgeList
