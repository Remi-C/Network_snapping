-----------------------------------------------------------
--
--Rémi-C , Thales IGN
--09/2014
--
--This script converts a topological network (edges beign polyline) to a topological network (edge being lines)
--
-----
--details
-----
---NOTE :  
--	doing so we want to conserve full history (what was the number of previous edges/node...)
--	
--warning : uses SFCGAL
-----------------------------------------------------------


CREATE SCHEMA IF NOT EXISTS network_for_snapping; 

SET search_path TO network_for_snapping, bdtopo_topological, bdtopo, topology, public; 



--break edges into pair of succesive points
--we give new node_id over 1000000 to nodes resulting from breaking

	DROP TABLE IF EXISTS new_successive_points; 
	CREATE TABLE new_successive_points AS 
	WITH input_data AS ( --getting all the edges we are going to break into segments .
		--no need to break edge if edge is already a segment and not a polyline!  (hence the tetst on point number)
		SELECT  * , ST_NumPoints(geom)  AS num_points
		 FROM edge_data
		 WHERE ST_NumPoints(geom) >2 
		--LIMIT  100 
	)
	 ,points AS ( --we break each edge into successiv points
		SELECT edge_id , dmp.path AS seg_id, dmp.geom ,num_points
		FROM input_data, ST_DumpPoints(geom) AS dmp 
	)
	 ,first_last AS ( --we get the information if the point is the first or the last of the edge. 
				--This is important, because if this is the case, the point will not be a new node but instead the old node that already exists betwen edges in bdtopo_topogical.node
		SELECT * , seg_id[1]= 1 AS is_first, seg_id[1]  = num_points AS is_last 
		FROM points
	) 
	,getting_old_node_id AS ( -- we genereate a new node id : if the point is already a node, keep the correct node_id (start_node or end_node).
			--else, generate a new id starting from 1000000
			-- the variable node_in_intersection keeps tracks of which point was already a node or not
		SELECT nni.edge_id, nni.seg_id[1] As seg_ordinality
			,   COALESCE(ed1.start_node, ed2.end_node, 1000000 + row_number() over()) as new_node_id
			, is_first OR is_last AS node_in_intersection
			,nni.geom
		FROM first_last as nni
			LEFT OUTER JOIN edge_data as ed1 ON (nni.edge_id = ed1.edge_id AND nni.is_first = TRUE)
			LEFT OUTER JOIN edge_data as ed2 ON (nni.edge_id = ed2.edge_id AND nni.is_last = TRUE)
	)
	SELECT *
	FROM getting_old_node_id ;


	DROP TABLE IF EXISTS nodes_for_output ; 
	CREATE TABLE nodes_for_output AS  
		--this is the sum of old nodes and newly created nodes. This should be output
			SELECT node_id
				,  ST_X( geom) AS X, ST_Y(geom) AS Y, ST_Z(geom) AS Z
				, 1 AS is_in_intersection
			FROM bdtopo_topological.node
		UNION ALL --in theory , we can use UNION ALL. This is a security to enforce no duplicates. We don't care about performance here
			SELECT new_node_id as node_id
				, ST_X( geom) AS X, ST_Y(geom) AS Y, ST_Z(geom) AS Z
				, 0 AS is_in_intersection
			FROM new_successive_points
			WHERE node_in_intersection = FALSE;
			--60 k lines
 

	DROP TABLE IF EXISTS edges_for_output ; 
	CREATE TABLE edges_for_output AS  
	WITH successive_points AS (--we get all the successive nodes from getting_old_node_id to create pair of successiv nodes
		SELECT edge_id, seg_ordinality 
			,  new_node_id as start_node
			, lead(new_node_id,1) OVER (PARTITION BY edge_id ORDER BY  seg_ordinality) AS end_node 
		FROM new_successive_points  	
	)
	,new_edges AS (--a pair of nodes is a new edge. We also get the width information.	
			-- we generate one new edge too much per old edge, we have to remove it using the "WHERE sp.end_node IS NOT NULL"
		SELECT 1000000+ row_number() over()  AS edge_id
			,sp.start_node
			,sp.end_node 
			,r.largeur AS width
			,sp.edge_id AS old_edge_id
		FROM successive_points AS sp
			LEFT OUTER JOIN edge_data AS ed ON (ed.edge_id = sp.edge_id)
			LEFT OUTER JOIN road AS  r ON (ed.ign_id = r.id)
		WHERE sp.end_node IS NOT NULL 
	)
		 SELECT ed.edge_id, ed.start_node, ed.end_node, r.largeur AS width
		FROM edge_data AS ed
		LEFT OUTER JOIN road AS  r ON (ed.ign_id = r.id)
	UNION --should be union all , security against duplicates.
		SELECT  edge_id
				, start_node
				, end_node 
				, width 
		FROM new_edges ;
	 


--now we restrain ourselve to an area for test :
	--this is delimitated by the table def_zone_export

	DROP TABLE IF EXISTS def_zone_export; 
	CREATE TABLE def_zone_export 
	( 	gid SERIAL PRIMARY KEY
		,id bigint
		,geom geometry(POLYGON,931008)
	);   
	INSERT INTO def_zone_export (id, geom) 
		VALUES( 1
			, ST_GeomFromText('POLYGON((650907.6  6860870.6 ,650956.9  6860895.8 ,651036.0  6860757.1 ,650983.7  6860736.6 ,650907.6  6860870.6 ))',931008) );


	DROP TABLE IF EXISTS nodes_for_output_in_export_area; 
	CREATE TABLE nodes_for_output_in_export_area AS 
		--we keep only the nodes inside the export area
	SELECT nou.*
	FROM nodes_for_output as nou, def_zone_export as dfz
	WHERE ST_WITHIN(ST_SetSRID(ST_MakePoint(nou.X,nou.Y,nou.Z),932011) , ST_Transform((dfz.geom),932011)) = TRUE;

	DROP TABLE IF EXISTS edges_for_output_in_export_area; 
	CREATE TABLE edges_for_output_in_export_area AS 
		--we keep only the edges whose both ends are exported nodes.
	SELECT eou.*
	FROM edges_for_output as eou , nodes_for_output_in_export_area AS nou1, nodes_for_output_in_export_area AS nou2
	WHERE eou.start_node = nou1.node_id AND eou.end_node = nou2.node_id  ;

	--importing the sidewalk observation (weighted points ) here
	DROP MATERIALIZED VIEW IF EXISTS weighted_sidewalk_observation_point ;
	CREATE MATERIALIZED VIEW weighted_sidewalk_observation_point AS 
	SELECT *
	FROM public.dblink( 
		 'hostaddr=127.0.0.1 port=5433 dbname=TerraMobilita user=postgres password=postgres'::text
		, 'SELECT * FROM cmm.lines_to_weighted_points '::text
		,TRUE)  AS   f(
			 qgis_id bigint,
			  area_id bigint,
			  seg_id integer,
			  geom geometry,
			  weight numeric
			 ) ; 

	--getting the observations inside the area 

	DROP TABLE IF EXISTS obs_for_output_in_export_area; 
	CREATE TABLE obs_for_output_in_export_area AS 
	WITH obs_in_area AS ( --these are the observations that are inside the defined area
			SELECT nou.*
			FROM weighted_sidewalk_observation_point AS  nou , def_zone_export as dfz
			WHERE ST_WITHIN(nou.geom,  dfz.geom ) = TRUE 
		)
	,edge_geom AS ( -- we reconstruct the edge geom to be able to assign observation to edges
		SELECT efo.*, ST_SetSRID(ST_MakeLine(ST_MakePoint(nfo1.X,nfo1.Y,nfo1.Z) ,ST_MakePoint(nfo2.X,nfo2.Y,nfo2.Z)  ),932011) as edge_geom
		FROM edges_for_output_in_export_area as efo
			LEFT JOIN nodes_for_output_in_export_area AS nfo1 ON (efo.start_node = nfo1.node_id)
			LEFT JOIN nodes_for_output_in_export_area AS nfo2 ON (efo.end_node = nfo2.node_id)
	)
	,map AS (--for each edge, we get observation closer than 5 meters
		SELECT DISTINCT ON (oia.qgis_id, eg.edge_id) oia.*, eg.edge_id
		FROM obs_in_area  AS oia ,edge_geom AS eg
		WHERE ST_DWithin( ST_Transform(oia.geom,932011), eg.edge_geom,5)=TRUE
		ORDER BY oia.qgis_id ASC, eg.edge_id ASC, ST_Distance(ST_Transform(oia.geom,932011),eg.edge_geom ) ASC
			
	)
	SELECT qgis_id AS obs_id , edge_id, ST_X(tgeom) AS X, ST_Y(tgeom) AS Y, COALESCE(ST_Z(tgeom),0) AS Z , 1::float AS confidence, COALESCE(weight,0) AS weight
	FROM map, ST_Transform( geom,932011) As tgeom ; 




	COPY  nodes_for_output_in_export_area TO '/media/sf_E_RemiCura/PROJETS/snapping/data/data_in_reduced_export_area/nodes_for_output_in_export_area.csv'
	 (FORMAT CSV, DELIMITER   ';') ;
	 COPY  edges_for_output_in_export_area TO '/media/sf_E_RemiCura/PROJETS/snapping/data/data_in_reduced_export_area/edges_for_output_in_export_area.csv'
	 (FORMAT CSV, DELIMITER   ';') ;
	 COPY  obs_for_output_in_export_area TO '/media/sf_E_RemiCura/PROJETS/snapping/data/data_in_reduced_export_area/obs_for_output_in_export_area.csv'
	 (FORMAT CSV, DELIMITER   ';') ;
	

	
	