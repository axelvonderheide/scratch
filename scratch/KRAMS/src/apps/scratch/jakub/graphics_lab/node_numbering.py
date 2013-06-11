       if point_numbers.any():
            for n, p in zip( point_numbers, points ):
                atext = tvtk.VectorText()
                atext.text=str(n)
                textMapper = tvtk.PolyDataMapper()
                textMapper.input_connection=atext.output_port

                textActor = tvtk.Follower(position=(p[0], p[1], p[2]), mapper=textMapper)
                textActor.scale=(0.2, 0.2, 0.2)
                textActor.property.color=(0, 0, 0)
                actors['number'+str(n)] = textActor
                #textActor.add_position=(0, -0.1, 0)