function [ClusterVertices] = ClusterBorder(ClusterMap)
%  CLUSTERBORDER function to draw line round a cluster of significant TFS
%  cells in colourmap e.g. using imagesc.
%  
% Takes a logical matrix showing a cluster of contiguous TFS cells, and
%  and finds the coordinates of the internal and external vertices of the
%  cluster, and orders them clockwise round the perimeter. The output of
%  ordered vertices can be used with patch.m to draw a line round the cluster
%% Pads the cluster array with zeros by one cell all round
Mask=abs(ClusterMap)>0;
PaddedMask=padarray(1*Mask,[1,1],0);
%% Finds linear and subscript indices of all cluster cells
LinIdx=find(PaddedMask==1);
[Row,Col]=ind2sub(size(PaddedMask),LinIdx);
ClusterSubs=[Row,Col];
%% Find neighbours of each cell in mask
NCellIdx=sub2ind(size(PaddedMask),Row-1,Col);
SCellIdx=sub2ind(size(PaddedMask),Row+1,Col);
WCellIdx=sub2ind(size(PaddedMask),Row,Col-1);
ECellIdx=sub2ind(size(PaddedMask),Row,Col+1);
NWCellIdx=sub2ind(size(PaddedMask),Row-1,Col-1);
NECellIdx=sub2ind(size(PaddedMask),Row-1,Col+1);
SWCellIdx=sub2ind(size(PaddedMask),Row+1,Col-1);
SECellIdx=sub2ind(size(PaddedMask),Row+1,Col+1);
ClusterCellNeighbours=[...
    PaddedMask(NWCellIdx),...
    PaddedMask(NCellIdx),...
    PaddedMask(NECellIdx),...
    PaddedMask(ECellIdx),...
    PaddedMask(SECellIdx),...
    PaddedMask(SCellIdx),...
    PaddedMask(SWCellIdx),...
    PaddedMask(WCellIdx)];
%% Find Edge cells
EdgeCellIdx=sum(ClusterCellNeighbours,2)<8; % Find all edge cells
EdgeSubs=ClusterSubs(EdgeCellIdx,:);
EdgeCellNeighbours=ClusterCellNeighbours(EdgeCellIdx,:);
ExtTL_Subs=EdgeSubs(sum(EdgeCellNeighbours(:,[2,8]),2)==0,:);
ExtTR_Subs=EdgeSubs(sum(EdgeCellNeighbours(:,[2,4]),2)==0,:);
ExtBL_Subs=EdgeSubs(sum(EdgeCellNeighbours(:,[6,8]),2)==0,:);
ExtBR_Subs=EdgeSubs(sum(EdgeCellNeighbours(:,[4,6]),2)==0,:);
IntTL_Subs=EdgeSubs(ismember(EdgeCellNeighbours(:,[8,1,2]), [1,0,1], 'rows'),:);
IntTR_Subs=EdgeSubs(ismember(EdgeCellNeighbours(:,2:4), [1,0,1], 'rows'),:);
IntBL_Subs=EdgeSubs(ismember(EdgeCellNeighbours(:,6:8), [1,0,1], 'rows'),:);
IntBR_Subs=EdgeSubs(ismember(EdgeCellNeighbours(:,4:6), [1,0,1], 'rows'),:);
%% Finds coordinates of all cluster corners
ExtTL_verts=[ExtTL_Subs(:,2)-0.5,ExtTL_Subs(:,1)-0.5];
ExtTR_verts=[ExtTR_Subs(:,2)+0.5,ExtTR_Subs(:,1)-0.5];
ExtBL_verts=[ExtBL_Subs(:,2)-0.5,ExtBL_Subs(:,1)+0.5];
ExtBR_verts=[ExtBR_Subs(:,2)+0.5,ExtBR_Subs(:,1)+0.5];
IntTL_verts=[IntTL_Subs(:,2)-0.5,IntTL_Subs(:,1)-0.5];
IntTR_verts=[IntTR_Subs(:,2)+0.5,IntTR_Subs(:,1)-0.5];
IntBL_verts=[IntBL_Subs(:,2)-0.5,IntBL_Subs(:,1)+0.5];
IntBR_verts=[IntBR_Subs(:,2)+0.5,IntBR_Subs(:,1)+0.5];
%% Makes column indicating of direction of next vertex (in polar coords)
ExtTLlabel=repmat(0,size(ExtTL_verts,1),1);
ExtTRlabel=repmat(-pi/2,size(ExtTR_verts,1),1);
ExtBLlabel=repmat(pi/2,size(ExtBL_verts,1),1);
ExtBRlabel=repmat(pi,size(ExtBR_verts,1),1);
IntTLlabel=repmat(pi/2,size(IntTL_verts,1),1);
IntTRlabel=repmat(0,size(IntTR_verts,1),1);
IntBLlabel=repmat(pi,size(IntBL_verts,1),1);
IntBRlabel=repmat(-pi/2,size(IntBR_verts,1),1);
Vertices_Cat=[...
    ExtTL_verts;...
    ExtTR_verts;...
    ExtBL_verts;...
    ExtBR_verts;...
    IntTL_verts;...
    IntTR_verts;...
    IntBL_verts;...
    IntBR_verts];
Directions_Cat=[...
    ExtTLlabel;...
    ExtTRlabel;...
    ExtBLlabel;...
    ExtBRlabel;...
    IntTLlabel;...
    IntTRlabel;...
    IntBLlabel;...
    IntBRlabel];
Vertices_Cat=[Vertices_Cat,Directions_Cat];
%% Order Vertices clockwise round perimeter of cluster
Nvertices=size(Vertices_Cat,1);
ThisVertex=Vertices_Cat(1,:);
Vertices_Cat(1,:)=[NaN,NaN,NaN];
StartVertex=ThisVertex;
PerimeterID=1;
iVertex=1;
ClusterVertices=NaN(Nvertices,3);
ClusterVertices(1,1:2)=ThisVertex(1:2);
ClusterVertices(iVertex,3)=PerimeterID;
for iVertex=2:Nvertices
    switch ThisVertex(3)
        case 0 % Go right
            NextIdx=find(Vertices_Cat(:,2)==ThisVertex(2));
            Possibles=Vertices_Cat(NextIdx,:);
            if StartVertex(2)==ThisVertex(2)
                Possibles=[StartVertex;Possibles];
            end
            Possibles=Possibles(Possibles(:,1)>ThisVertex(1),:);
            NextVertex=Possibles(Possibles(:,1)==min(Possibles(:,1)),:);
            NextVertex=NextVertex(1,:); % to handle cases where the same vertex occurs twice
        case pi % Go left
            NextIdx=find(Vertices_Cat(:,2)==ThisVertex(2));
            Possibles=Vertices_Cat(NextIdx,:);
            if StartVertex(2)==ThisVertex(2)
                Possibles=[StartVertex;Possibles];
            end
            Possibles=Possibles(Possibles(:,1)<ThisVertex(1),:);
            NextVertex=Possibles(Possibles(:,1)==max(Possibles(:,1)),:);
            NextVertex=NextVertex(1,:);
        case pi/2 % Go up
            NextIdx=find(Vertices_Cat(:,1)==ThisVertex(1));
            Possibles=Vertices_Cat(NextIdx,:);
            if StartVertex(1)==ThisVertex(1)
                Possibles=[StartVertex;Possibles];
            end
            Possibles=Possibles(Possibles(:,2)<ThisVertex(2),:);
            NextVertex=Possibles(Possibles(:,2)==max(Possibles(:,2)),:);
            NextVertex=NextVertex(1,:);
        case -pi/2 % Go down
            NextIdx=find(Vertices_Cat(:,1)==ThisVertex(1));
            Possibles=Vertices_Cat(NextIdx,:);
            if StartVertex(1)==ThisVertex(1)
                Possibles=[StartVertex;Possibles];
            end
            Possibles=Possibles(Possibles(:,2)>ThisVertex(2),:);
            NextVertex=Possibles(Possibles(:,2)==min(Possibles(:,2)),:);
            NextVertex=NextVertex(1,:);
    end
    if ~ismember(NextVertex,StartVertex,"rows")
        ThisVertex=NextVertex;
    else
        PerimeterID=PerimeterID+1;
        VerticesLeftIdx=find(~isnan(Vertices_Cat(:,1)));
        ThisVertex=Vertices_Cat(VerticesLeftIdx(1),:);
        StartVertex=ThisVertex;
    end
    [~,idx]=ismember(ThisVertex,Vertices_Cat,'rows');
    Vertices_Cat(idx,:)=[NaN,NaN,NaN];
    ClusterVertices(iVertex,1:2)=ThisVertex(1:2);
    ClusterVertices(iVertex,3)=PerimeterID;
end
%% Remove padding from vertex coordinates
ClusterVertices(:,1:2)=ClusterVertices(:,1:2)-1;
end