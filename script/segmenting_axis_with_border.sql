

SET search_path TO network_for_snapping, bdtopo_topological, bdtopo, topology,rc_lib, public; 

	DROP TABLE IF EXISTS test_splitting_edge ; 
	CREATE TABLE  test_splitting_edge AS
	WITh edges_to_split AS (
		SELECT -- count(*) -- edge_id, intersection_limit1, intersection_limit2
			edge_id, ed.geom::geometry(linestringZ,932011),  intersection_limit1::geometry(point,932011), intersection_limit2::geometry(point,932011)
		FROM street_amp.result_axis AS ra 
			LEFT OUTER JOIN bdtopo_topological.edge_data AS ed USING (edge_id)
			LEFT OUTER JOIN bdtopo_topological.node AS n1 ON (n1.node_id = ed.start_node)
			LEFT OUTER JOIN bdtopo_topological.node AS n2 ON (n2.node_id = ed.end_node)
			--, (SELECT min(ST_Distance(dmp.geom, )
		WHERE intersection_limit1 IS NOT NULL AND intersection_limit2 IS NOT NULL --25071
			AND ST_NumPoints(ed.geom) = 2
			AND ST_Length(ra.section2) > 9
			AND in_n_intersection[2] = 0 
			AND ST_DWIthin(n1.geom,ra.intersection_limit1,5)= FALSE 
			AND ST_DWIthin(n2.geom,ra.intersection_limit1,5)= FALSE 
			AND ST_DWIthin(n1.geom,ra.intersection_limit2,5)= FALSE 
			AND ST_DWIthin(n2.geom,ra.intersection_limit2,5)= FALSE 
	)
	, splitting AS (
		SELECT edge_id, ST_AddPoint(geom,npoint,1)::geometry(linestringZ,932011) AS new_geom
		FROM edges_to_split
			, ST_LineLocatePoint(geom, intersection_limit1) AS curvabs
			, ST_LineInterpolatePoint(geom, curvabs) AS npoint  
	)
	UPDATE  bdtopo_topological.edge_data AS ed  SET geom=
	 s.new_geom
	FROM splitting AS s
	WHERE s.edge_id = ed.edge_id
 
LIMIT 1

SELECT  count(*) -- edge_id, intersection_limit1, intersection_limit2
FROM street_amp.result_axis AS ra 
LEFT OUTER JOIN bdtopo_topological.edge_data AS ed USING (edge_id)
WHERE ST_NumPoints(ed.geom) = 2