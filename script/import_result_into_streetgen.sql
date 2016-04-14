

CREATE EXTENSION IF NOT EXISTS dblink; 
CREATE SCHEMA IF NOT EXISTS network_for_snapping; 
SET search_path TO network_for_snapping, bdtopo_topological, bdtopo, topology,rc_lib, public; 



--importing data
DROP MATERIALIZED VIEW IF EXISTS optimized_edges ;
	CREATE MATERIALIZED VIEW optimized_edges AS 
	SELECT edge_id, start_node, end_node, width, old_width, cost_edge, iteration, old_edge_id, n_obs, old_geom::geometry(linestringZ,932011), geom::geometry(linestringZ,932011)
	FROM dblink( 
		 'hostaddr=127.0.0.1 port=5433 dbname=street_gen_3 user=postgres password=postgres'::text
		, '
		WITH obs AS (
			SELECT edge_id,  count(*) n_obs
			FROM network_for_snapping.obs_for_output_in_export_area AS ou
			GROUP BY edge_id
		)
		SELECT eo.edge_id, start_node, end_node, oe.width, eo.width AS old_width, cost_edge, iteration, eo.old_edge_id, n_obs , eo.geom AS old_geom, oe.geom 
		FROM network_for_snapping.optimized_edges AS oe
			LEFT OUTER JOIN network_for_snapping.edges_for_output as eo ON (eo.edge_id = oe.gid)  
			LEFT OUTER JOIN obs ON (obs.edge_id = oe.gid) ; '::text
		,TRUE)  AS   f(
			 edge_id bigint,
			  start_node bigint,
			  end_node bigint,
			  width float,
			  old_width float,
			  cost_edge float,
			  iteration int,
			  old_edge_id int, 
			  n_obs int, 
			  old_geom geometry(linestringZ,932011) ,
			  geom geometry(linestringZ,932011) 
			 ) ; 
CREATE INDEX ON optimized_edges (edge_id) ;
CREATE INDEX ON optimized_edges (old_edge_id) ;

--reconstructing node table:

DROP TABLE IF EXISTS rec_node ;
CREATE TABLE rec_node  AS 
SELECT node_id, (array_agg(node_geom))[1]::geometry(pointz,932011) AS node_geom, sum(n_obs) as n_obs
FROM (
	SELECT start_node AS node_id, ST_StartPoint(geom) as node_geom, n_obs
	FROM optimized_edges
	UNION ALL
	SELECT end_node AS node_id, ST_EndPoint(geom) as node_geom , n_obs
	FROM optimized_edges )
	AS sub
GROUP BY node_id;
ALTER TABLE rec_node ADD PRIMARY KEY(node_id) ;


DROP TABLE IF EXISTS rec_node_old ;
CREATE TABLE rec_node_old  AS 
SELECT node_id, (array_agg(old_node_geom))[1]::geometry(pointz,932011) As old_node_geom, sum(n_obs) as n_obs
FROM (
	SELECT start_node AS node_id, ST_StartPoint(old_geom) as old_node_geom, n_obs
	FROM optimized_edges
	UNION ALL
	SELECT end_node AS node_id, ST_EndPoint(old_geom) as old_node_geom , n_obs
	FROM optimized_edges ) AS sub
	GROUP BY node_id;
ALTER TABLE rec_node_old ADD PRIMARY KEY(node_id) ; 

 

----- how much change was applied on network? 
SELECT avg_node_change, avg_width_change
FROM (
SELECT avg(COALESCE(ST_Distance(old_node_geom,  node_geom),0)) AS avg_node_change
FROM rec_node
	LEFT OUTER JOIN rec_node_old USING(node_id) ) AS sub
, ( SELECT avg(abs(old_width-width)) avg_width_change
FROM optimized_edges 
) as sub2  ; 
--1.31314996900513	1.42093499240389

-- regrouping by old edge_id, for each old_edge_id, we got a successive list oe edge_id, each with n_obs, and a width
-- goal is to regroup edge with a compatible width
DROP TABLE IF EXISTS edge_with_label ;
CREATE TABLE edge_with_label AS
	WITH edges_for_one_old_edge AS (
		SELECT edge_id, old_edge_id, width, n_obs
		FROM optimized_edges
		WHERE n_obs IS NOT NULL
			AND n_obs >= 4
		--AND old_edge_id   <100
		ORDER BY edge_id
	)
	 , fabricating_data AS (
		SELECT old_edge_id, array_agg(edge_id::int ORDER BY edge_id) AS edge_ids, array_agg(width ORDER BY edge_id) AS widths,  array_agg(n_obs ORDER BY edge_id) AS n_obss
		FROM edges_for_one_old_edge 
		GROUP BY old_edge_id 
	)
	, labels_base AS(
	SELECT f.gid AS edge_id, f.label , old_edge_id
	FROM fabricating_data
		, rc_lib.rc_py_cluster_points_dbscan( edge_ids, widths  -- , n_obss
			,ndim:=1,  max_dist := 0.6, min_sample:=1) AS f
	)
	SELECT edge_id, label, lb.old_edge_id, width, geom::geometry(linestringz,932011)
	FROM labels_base AS lb
		LEFT OUTER JOIN optimized_edges AS oe USING (edge_id);

 
----- for edges of one old_edge_id, we have an ordered list of edge id and an ordered list of label. 
-- non adjacent edges should not have the same label
-- to enforce this, we use the trick of row_number() over(partition by old_edge_id, label order by edge_id)
-- the new table contains for each edge with sufficient observation its new label
	DROP TABLE IF EXISTS n_label ;
	CREATE TABLE n_label  AS   
	WITH i_data AS (
		SELECT *
		FROM edge_with_label
		--WHERE old_edge_id = 470
		ORDER BY old_edge_id 
	)
	 , new_label AS (
	SELECT old_edge_id, edge_id, label,width,  first_value(edge_id) over(partition by old_edge_id, label, rel_ordering ORDER BY old_edge_id, edge_id)  AS n_label
	FROM (
		SELECT old_edge_id, edge_id, label, width, row_number() over(ORDER BY old_edge_id, edge_id) - row_number() over(partition by old_edge_id, label order by old_edge_id, edge_id) AS rel_ordering
		FROM i_data
	) AS sub
	)
	SELECT edge_id, old_edge_id, n_label
	FROM new_label ;
	CREATE INDEX ON n_label (edge_id); 
	CREATE INDEX ON n_label (old_edge_id);  

-- regrouping edges with same label, computing weighted median of width, geometry merging of same label

DROP TABLE IF EXISTS line_sewing ;
CREATE TABLE line_sewing  AS 
SELECT n_label,  nl.old_edge_id,  ST_LineMerge(ST_Collect(geom order by edge_id)) as n_edge_geom-- ::geometry(multilinestringZ,932011)
	--, array_agg(edge_id order by edge_id) as edge_ids
	,  rc_py_weighted_median( array_agg(width order by edge_id), array_agg(n_obs order by edge_id)  ) as wmed 
	, sum(n_obs) AS n_obs
FROM n_label AS nl
	LEFT OUTER JOIN optimized_edges AS oe USING (edge_id)
GROUP BY nl.old_edge_id, n_label ; 

ALTER TABLE line_sewing ALTER COLUMN n_edge_geom  TYPE  geometry( linestringZ,932011) USING ST_GeometryN(n_edge_geom,1) ;
CREATE INDEX ON line_sewing USING GIST(n_edge_geom) ; 

SELECT *
FROM  line_sewing
WHERE ST_IsVAlid(n_edge_geom) = FALSE OR ST_IsSimple(n_edge_geom) = FALSE ; 

--now dealing with edges with too few observation to be judged.
-- their widht is assigned to the closest other edge width witht the max observation, and is old_width else

DROP TABLE IF EXISTS no_obs_edge_with_n_width ; 
CREATE TABLE no_obs_edge_with_n_width AS 
-- WITH edge_not_labelled_with_width AS (
	(SELECT DISTINCT ON (oe.old_edge_id, edge_id )  
		oe.edge_id,  oe.old_edge_id , ls.n_label AS new_edge_id_map
		, COALESCE(ls.wmed, oe.width) AS n_width
		, oe.geom::geometry(linestringZ,932011) 
	FROM optimized_edges AS oe
		LEFT OUTER JOIN n_label AS nl USING (edge_id )
		, line_sewing AS ls
	WHERE nl.n_label IS NULL
		AND oe.old_edge_id = ls.old_edge_id 
		-- AND ST_DWithin(ls.n_edge_geom, nl.geom,100)=true
	ORDER BY old_edge_id, edge_id,  ST_Distance(ls.n_edge_geom, oe.geom) , ls.n_obs DESC  )
	UNION ALL
	SELECT oe.edge_id,  oe.old_edge_id , oe.edge_id AS new_edge_id_map 
		, oe.width AS n_width
		, oe.geom::geometry(linestringZ,932011) 
	FROM optimized_edges AS oe
		LEFT OUTER JOIN n_label AS nl USING (edge_id )
	WHERE NOT EXISTS (
		SELECT 1 
		FROM line_sewing AS ls
		WHERE  ls.old_edge_id = oe.old_edge_id
	) ; 

---- now sewing no observation edges with other edges 
DROP TABLE IF EXISTS all_edges_sewed ; 
CREATE TABLE all_edges_sewed AS
	SELECT old_edge_id, n_label AS new_edge_id
		,  ST_LineMerge(
			ST_Collect(n_edge_geom  order by edge_id)
			)  AS n_geom
		, (array_agg(wmed))[1] AS wmed, sum(n_obs) AS n_obs
	FROM (
	SELECT  old_edge_id, n_label,n_label AS edge_id, n_edge_geom, wmed, n_obs
	FROM line_sewing AS ls
	UNION ALL 
	SELECT old_edge_id, new_edge_id_map AS n_label, edge_id , geom,n_width AS wmed, 0 AS n_obs
	FROM no_obs_edge_with_n_width AS noe
	--LIMIT 100
	) AS sub
	GROUP BY old_edge_id, n_label ; 

ALTER TABLE all_edges_sewed ALTER COLUMN n_geom  TYPE  geometry( linestringZ,932011) USING ST_MakeValid(ST_GeometryN(n_geom,1)) ;
CREATE INDEX ON all_edges_sewed USING GIST(n_geom) ; 
ALTER TABLE all_edges_sewed ADD PRIMARY KEY (new_edge_id) ; 
DELETE FROM all_edges_sewed WHERE ST_IsEmpty(n_geom)=TRUE OR ST_IsValid(n_geom) = FALSE OR ST_IsSImple(n_geom) = FALSE
	OR wmed IS NULL OR wmed <=0 ; 

-- DROP TABLE IF EXISTS temp_visu_invalid ; 
-- CREATE TABLE temp_visu_invalid AS 
-- SELECT n_label, old_edge_id, n_geom::geometry( linestringZ,932011), wmed, n_obs
-- FROM all_edges_sewed 
-- WHERE ST_IsEmpty(n_geom)=TRUE OR ST_IsValid(n_geom) = FALSE OR ST_IsSImple(n_geom) = FALSE
-- 	OR wmed IS NULL OR wmed <=0 ;



-- now we have a linestring layer repersenting the road axis, ready to be topologised
-- yet, we still lack other initial data, therefore we are going ot map the new edges to the old ones, so we are able to transfer attributes such as lane number and so .

-- importing BDTOPO.route, then converting it into bdtopo roadwith streetgen
--ALTER TABLE bdtopo.road RENAME TO road_ori ; 
--mapping between all_edges_sewed and road


	DROP TABLE IF EXISTS final_optimized_edges ; 
	CREATE TABLE final_optimized_edges AS 
	SELECT DISTINCT ON (ae.new_edge_id) 
	row_number() over() as gid, new_edge_id::text AS id, 
		nature  ,
		nom_voie_g   
		,nom_voie_d 
		,importance ,--required
		fictif ,--required
		franchisst , 
		wmed as largeur ,--required
			nb_voies  ,--required
			pos_sol   ,
			sens  ,
			etat  ,
			z_ini  ,
			z_fin  ,
			NULL AS radius ,
			ST_Multi(ST_transform(n_geom,931008) )AS geom
	FROM all_edges_sewed AS ae , bdtopo.route AS ro 
	WHERE ST_DWIthin(ae.n_geom, St_Transform(ro.geom,932011),20) = TRUE
	ORDER BY new_edge_id, ST_Distance(ae.n_geom, St_Transform(ro.geom,932011)), ST_Area(ST_INtersection(ST_Buffer(ae.n_geom,3),ST_Buffer(St_Transform(ro.geom,932011), 3))) DESC, ro.id ; 

	CREATE INDEX ON final_optimized_edges USING GIST(geom) ; 
--now remove edges that may intersects others (eliminating common nodes) 
	DROP TABLE IF EXISTS  error_in_network ;
	CREATE TABLE error_in_network AS
	SELECT row_number()over() as gid
		, case when ST_Length(e1.geom) > st_length(e2.geom) then e2.gid else e1.gid END as gid_to_be_deleted
		, cent::geometry(point,931008)
		, (CASE WHEN ST_Length(e1.geom) > st_length(e2.geom) THEN e2.geom ELSE e1.geom END)::geometry(multilinestringZ,931008) AS inital_geom
	FROM final_optimized_edges AS e1
		, final_optimized_edges AS e2 
		,ST_Centroid(ST_Intersection(e1.geom,e2.geom)) as cent
	WHERE e1.gid < e2.gid 
		AND ST_Intersects(e1.geom,e2.geom) = TRUE 
		AND ( abs(ST_LineLocatePoint(ST_GeometryN(e1.geom,1),cent)*2-1) BETWEEN 0 AND 0.9
		OR abs(ST_LineLocatePoint( ST_GeometryN(e2.geom,1),cent)*2-1) BETWEEN 0 AND 0.9) ;

/*	WITH to_be_deleted AS (
		SELECT distinct gid_to_be_deleted
		FROM error_in_network
	)
	DELETE FROM final_optimized_edges AS fo USING  to_be_deleted AS tbd 
	WHERE fo.gid =  tbd.gid_to_be_deleted ; 
*/	


	-- ALTER TABLE bdtopo.route ALTER COLUMN geom type  geometry(multilinestringZm,931008) using st_setsrid(geom,931008) ; 
	--ALTER TABLE final_optimized_edges rename to road; 
	--ALTER TABLE road SET SCHEMA bdtopo ; 
	--ALTER TABLE road ADD PRIMARY KEY(gid) ; 
	-- create index ! see streetgen 


    --now : finish isntall and generate streetgen